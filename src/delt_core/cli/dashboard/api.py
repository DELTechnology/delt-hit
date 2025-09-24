from pathlib import Path
import pandas as pd
import dash
from dash import dcc, html, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import yaml
import re
from delt_core.utils import read_yaml


config_path = Path(
    '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml')
cfg = read_yaml(config_path)
counts_path = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/selections/AG24_1/counts.txt')

def load_config(config_path):
    try:
        with open(config_path, 'r') as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        return {
            'experiment': {
                'name': 'test-1',
                'fastq_path': '~/data/DECLT-DB/fastq_files/368061_1-241105_AG_BZ_NC_pool1_NF_S3_R1_001.fastq.gz',
                'save_dir': '~/data/DECLT-DB/experiments',
                'num_cores': 10
            },
            'selections': {
                'AG24_1': {
                    'operator': 'A. Gloger',
                    'date': '2024-09-26',
                    'target': 'No Protein',
                    'group': 'no_protein',
                    'beads': 'Dynabeads SA C1',
                    'blocking': 'Biotin',
                    'buffer': 'PBS-T',
                    'protocol': 'DECL_5W',
                    'S0': 'ACACAC',
                    'S1': 'CGCTCGATA'
                }
            }
        }

def load_counts(counts_path):
    try:
        return pd.read_csv(counts_path, sep='\t')
    except FileNotFoundError:
        return pd.DataFrame({
            'code_1': [1, 1, 1, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 20],
            'code_2': [1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 20],
            'count':  [6, 3, 14, 15, 6, 6, 4, 15, 3, 2, 4, 8, 6, 3]
        })

def get_available_codes(df):
    return [col for col in df.columns if col.startswith('code_')]

def marginalize_counts(df, selected_codes):
    if not selected_codes:
        return df
    codes_to_use = selected_codes[:3] if len(selected_codes) > 3 else selected_codes
    grouped = df.groupby(codes_to_use, dropna=False)['count'].sum().reset_index()
    return grouped

def parse_code_ranges(range_str: str, code_cols: list[str]) -> dict:
    """
    Parse a string like '1-3;2-5,8-10;3-5' into a dict:
    {'code_1': {1,2,3}, 'code_2': {2,3,4,5,8,9,10}, 'code_3': {3,4,5}}
    Extra segments are ignored; missing segments mean 'no filtering' for that code.
    """
    if not range_str or not isinstance(range_str, str):
        return {}

    segments = [seg.strip() for seg in range_str.split(';')]
    filters = {}

    for idx, seg in enumerate(segments):
        if idx >= len(code_cols):
            break
        if not seg:
            continue
        allowed = set()
        for part in seg.split(','):
            part = part.strip()
            if not part:
                continue
            m = re.match(r'^(-?\d+)\s*-\s*(-?\d+)$', part)
            if m:
                a, b = int(m.group(1)), int(m.group(2))
                lo, hi = (a, b) if a <= b else (b, a)
                allowed.update(range(lo, hi + 1))
            else:
                # single integer
                try:
                    allowed.add(int(part))
                except ValueError:
                    # ignore malformed tokens
                    pass
        if allowed:
            filters[code_cols[idx]] = allowed
    return filters

def apply_code_filters(df: pd.DataFrame, filters: dict) -> pd.DataFrame:
    """Keep rows where each filtered code column is in the allowed set."""
    if not filters:
        return df
    mask = pd.Series(True, index=df.index)
    for col, allowed in filters.items():
        if col in df.columns:
            mask &= df[col].isin(allowed)
    return df.loc[mask]

def create_config_cards(config):
    cards = []
    exp_info = config.get('experiment', {})
    exp_items = [
        html.P([html.Strong("Name: "), exp_info.get('name', 'N/A')]),
        html.P([html.Strong("FASTQ Path: "), html.Code(exp_info.get('fastq_path', 'N/A'))]),
        html.P([html.Strong("Save Directory: "), html.Code(exp_info.get('save_dir', 'N/A'))]),
        html.P([html.Strong("CPU Cores: "), str(exp_info.get('num_cores', 'N/A'))])
    ]
    cards.append(
        html.Div([
            html.H3("Experiment Configuration", className="text-lg font-bold mb-2"),
            html.Div(exp_items)
        ], className="bg-white p-4 rounded-lg shadow mb-4")
    )
    selections = config.get('selections', {})
    for selection_id, selection_info in selections.items():
        selection_items = [html.P([html.Strong(f"{k.replace('_', ' ').title()}: "), str(v)])
                           for k, v in selection_info.items()]
        cards.append(
            html.Div([
                html.H3(f"Selection: {selection_id}", className="text-lg font-bold mb-2"),
                html.Div(selection_items)
            ], className="bg-white p-4 rounded-lg shadow mb-4")
        )
    return cards


def dashboard(*, config_path: Path, counts_path: Path):

    # Load data
    config = read_yaml(config_path)
    counts_df = pd.read_csv(counts_path, sep='\t')

    selection_name = counts_path.parent.name
    (_, selection), = list(filter(lambda x: x[0] == selection_name, config['selections'].items()))
    config['selections'] = {selection_name: selection}

    available_codes = get_available_codes(counts_df)

    # Defaults for filters
    def default_code_range_string(df, code_cols):
        segs = []
        for c in code_cols:
            if c in df and df[c].notna().any():
                try:
                    mn, mx = int(df[c].min()), int(df[c].max())
                    segs.append(f"{mn}-{mx}")
                except Exception:
                    segs.append("")  # non-integer codes: leave blank
            else:
                segs.append("")
        return ';'.join(segs) if segs else ""

    # Initialize Dash app
    app = dash.Dash(__name__)

    axis_options = [{'label': 'None', 'value': 'None'}] + \
                   [{'label': c, 'value': c} for c in available_codes] + \
                   [{'label': 'count', 'value': 'count'}]

    app.layout = html.Div([
        html.Div([
            html.H1("DECL Experiment Dashboard",
                    className="text-3xl font-bold text-center mb-8 text-blue-600")
        ]),

        html.Div([
            html.H2("Configuration", className="text-2xl font-bold mb-4"),
            html.Div(create_config_cards(config), id="config-cards")
        ], className="mb-8"),

        # Controls
        html.Div([
            html.H2("Counts Analysis", className="text-2xl font-bold mb-4"),
            html.Div([
                html.Div([
                    html.Label("X-axis:", className="font-semibold mb-2 block"),
                    dcc.Dropdown(
                        id='x-axis-selector',
                        options=axis_options,
                        value=available_codes[0] if available_codes else 'None',
                        placeholder="Select X-axis variable..."
                    )
                ], className="w-1/4 pr-2"),

                html.Div([
                    html.Label("Y-axis:", className="font-semibold mb-2 block"),
                    dcc.Dropdown(
                        id='y-axis-selector',
                        options=axis_options,
                        value='count',
                        placeholder="Select Y-axis variable..."
                    )
                ], className="w-1/4 px-2"),

                html.Div([
                    html.Label("Z-axis (3D):", className="font-semibold mb-2 block"),
                    dcc.Dropdown(
                        id='z-axis-selector',
                        options=axis_options,
                        value='None',
                        placeholder="Select Z-axis variable (select to enable 3D)..."
                    )
                ], className="w-1/4 px-2"),

                html.Div([
                    html.Label("Options:", className="font-semibold mb-2 block"),
                    dcc.Checklist(
                        id='color-by-count',
                        options=[{'label': 'Color by count', 'value': 'on'}],
                        value=[],
                        inputStyle={"margin-right": "6px"}
                    ),
                    dcc.Checklist(
                        id='size-by-count',
                        options=[{'label': 'Size by count', 'value': 'on'}],
                        value=[],
                        inputStyle={"margin-right": "6px"},
                        style={"margin-top": "8px"}
                    ),
                ], className="w-1/4 pl-2"),
            ], className="flex bg-white p-4 rounded-lg shadow mb-4")
        ]),

        # Filters
        html.Div([
            html.H2("Filters", className="text-2xl font-bold mb-4"),
            html.Div([
                html.Div([
                    html.Label("Code Ranges (e.g. '1-3;2-5,8-10;3-5')", className="font-semibold mb-2 block"),
                    # Code Ranges
                    dcc.Input(
                        id='filter-codes',
                        type='text',
                        value=default_code_range_string(counts_df, available_codes),
                        style={'width': '100%'},
                        className="mb-2 border border-gray-300 rounded-md px-2 py-1 "
                                  "focus:outline-none focus:ring-1 focus:ring-blue-600"
                    ),
                ], className="w-1/2 pr-2"),

                html.Div([
                    html.Label("Min Count", className="font-semibold mb-2 block"),
                    # Min Count
                    dcc.Input(
                        id='filter-min-count',
                        type='number',
                        value=int(counts_df['count'].min()) if not counts_df.empty else 0,
                        step=1,
                        style={'width': '100%'},
                        className="mb-2 border border-gray-300 rounded-md px-2 py-1 "
                                  "focus:outline-none focus:ring-1 focus:ring-blue-600"
                    ),
                ], className="w-1/4 px-2"),

                html.Div([
                    html.Label("Max Count", className="font-semibold mb-2 block"),
                    # Max Count
                    dcc.Input(
                        id='filter-max-count',
                        type='number',
                        value=int(counts_df['count'].max()) if not counts_df.empty else 0,
                        step=1,
                        style={'width': '100%'},
                        className="mb-2 border border-gray-300 rounded-md px-2 py-1 "
                                  "focus:outline-none focus:ring-1 focus:ring-blue-600"
                    ),
                ], className="w-1/4 pl-2"),
            ], className="flex bg-white p-4 rounded-lg shadow mb-2"),

            html.Div([
                html.Button("Filter", id='filter-button', n_clicks=0,
                            className="px-4 py-2 bg-blue-600 text-white rounded shadow"),
                html.Button("Reset", id='reset-button', n_clicks=0,
                            className="px-4 py-2 bg-gray-200 text-gray-800 rounded shadow ml-2"),
            ], className="mb-4")
        ], className="mb-8"),

        # Stats + Plot
        html.Div([
            html.H3("Dataset Statistics", className="text-xl font-bold mb-2"),
            html.Div(id="stats-display", className="bg-gray-100 p-4 rounded-lg mb-4")
        ]),
        html.Div([
            dcc.Graph(id='counts-plot', style={'height': '800px', 'width': '100%'})
        ], className="bg-white p-4 rounded-lg shadow")

    ], className="container mx-auto p-4 bg-gray-50 min-h-screen")

    # Reset callback: restore defaults into the three filter inputs
    @app.callback(
        [Output('filter-codes', 'value'),
         Output('filter-min-count', 'value'),
         Output('filter-max-count', 'value')],
        [Input('reset-button', 'n_clicks')],
        prevent_initial_call=True
    )
    def reset_filters(n_clicks):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        return (
            default_code_range_string(counts_df, available_codes),
            int(counts_df['count'].min()) if not counts_df.empty else 0,
            int(counts_df['count'].max()) if not counts_df.empty else 0
        )

    @app.callback(
        [Output('counts-plot', 'figure'),
         Output('stats-display', 'children')],
        [Input('x-axis-selector', 'value'),
         Input('y-axis-selector', 'value'),
         Input('z-axis-selector', 'value'),
         Input('color-by-count', 'value'),
         Input('size-by-count', 'value'),
         Input('filter-button', 'n_clicks')],
        [State('filter-codes', 'value'),
         State('filter-min-count', 'value'),
         State('filter-max-count', 'value')],
        prevent_initial_call=False
    )
    def update_plot_and_stats(x_axis, y_axis, z_axis, color_by_count, size_by_count,
                              n_clicks_filter, range_str, min_count, max_count):
        # Require at least one axis
        if (x_axis == 'None') and (y_axis == 'None') and (z_axis == 'None'):
            empty_fig = go.Figure()
            empty_fig.add_annotation(text="Please select at least one axis",
                                     xref="paper", yref="paper", x=0.5, y=0.5,
                                     showarrow=False, font_size=16)
            return empty_fig, "No axes selected"

        # 1) Apply code range filtering on the raw data
        code_filters = parse_code_ranges(range_str, available_codes)
        df_filtered = apply_code_filters(counts_df, code_filters)

        # 2) Collect selected code columns (from axes) to aggregate
        selected_codes = []
        for axis in (x_axis, y_axis, z_axis):
            if axis not in ('None', 'count') and axis is not None:
                selected_codes.append(axis)
        selected_codes = list(dict.fromkeys(selected_codes))

        # 3) Aggregate/marginalize
        if selected_codes:
            plot_data = marginalize_counts(df_filtered, selected_codes)
        else:
            plot_data = pd.DataFrame({'count': [df_filtered['count'].sum()]})

        # 4) Apply min/max count filters on aggregated data (if provided)
        if min_count is not None:
            plot_data = plot_data.loc[plot_data['count'] >= int(min_count)]
        if max_count is not None:
            plot_data = plot_data.loc[plot_data['count'] <= int(max_count)]

        # Compute stats
        total_entries = len(plot_data)
        total_counts = int(plot_data['count'].sum()) if 'count' in plot_data else 0
        avg_count = float(plot_data['count'].mean()) if 'count' in plot_data and total_entries else 0.0
        max_count_val = int(plot_data['count'].max()) if 'count' in plot_data and total_entries else 0
        min_count_val = int(plot_data['count'].min()) if 'count' in plot_data and total_entries else 0
        codes_used_text = f"Codes used (grouped): {', '.join(selected_codes)}" if selected_codes else "No codes (total aggregation)"

        stats = html.Div([
            html.P(codes_used_text),
            html.P(f"Data points (after filters): {total_entries}"),
            html.P(f"Total counts (after filters): {total_counts:,}"),
            html.P(f"Average count: {avg_count:.1f}"),
            html.P(f"Count range: {min_count_val} - {max_count_val}"),
            html.P(f"Applied code filters: {range_str or '(none)'}"),
            html.P(
                f"Min/Max count filters: {min_count if min_count is not None else '-'} / {max_count if max_count is not None else '-'}")
        ])

        # Visual encodings
        color_dim = 'count' if ('on' in (color_by_count or [])) else None
        size_dim = 'count' if ('on' in (size_by_count or [])) else None

        use_3d = (z_axis != 'None')

        def labels_for(x=None, y=None, z=None):
            lab = {}
            if x: lab[x] = x
            if y: lab[y] = y
            if z: lab[z] = z
            lab['count'] = 'Count'
            return lab

        fig = go.Figure()

        # 3D plot
        if use_3d:
            if x_axis == 'None' or y_axis == 'None':
                fig.add_annotation(text="For 3D, please select X, Y, and Z.",
                                   xref="paper", yref="paper", x=0.5, y=0.5,
                                   showarrow=False, font_size=16)
            else:
                fig = px.scatter_3d(
                    plot_data,
                    x=x_axis if x_axis != 'None' else None,
                    y=y_axis if y_axis != 'None' else None,
                    z=z_axis if z_axis != 'None' else None,
                    color=color_dim,
                    size=size_dim,
                    labels=labels_for(x_axis, y_axis, z_axis),
                    title=f"3D: {x_axis} vs {y_axis} vs {z_axis}"
                          f"{' • color=count' if color_dim else ''}"
                          f"{' • size=count' if size_dim else ''}"
                )

                # Make 3D fill space better and avoid a squashed z axis
                fig.update_layout(
                    scene=dict(
                        aspectmode='cube',  # equal on-screen scale for x/y/z
                        xaxis=dict(zeroline=False),
                        yaxis=dict(zeroline=False),
                        zaxis=dict(zeroline=False),
                        domain=dict(x=[0.0, 1.0], y=[0.0, 1.0])  # occupy full area
                    )
                )

                # Make sized points more visible in 3D
                # fig.update_traces(
                #     marker=dict(
                #         sizemode='diameter',
                #         sizeref=None,  # let Plotly auto-compute
                #     ),
                #     selector=dict(type='scatter3d')
                # )

                # If not sizing by count, set a comfortable fixed size in 3D
                if size_dim is None:
                    fig.update_traces(marker=dict(size=3), selector=dict(type='scatter3d'))
        else:
            # 2D cases
            if x_axis == 'None' or y_axis == 'None':
                if x_axis != 'None':
                    if x_axis == 'count':
                        fig = px.histogram(
                            plot_data, x='count',
                            title='Distribution of Counts',
                            labels={'count': 'Count', 'count_count': 'Frequency'}
                        )
                    else:
                        fig = px.bar(
                            plot_data, x=x_axis, y='count',
                            title=f'Counts by {x_axis}',
                            labels={'count': 'Count', x_axis: x_axis}
                        )
                elif y_axis != 'None':
                    if y_axis == 'count':
                        fig = px.bar(x=['Total'], y=[total_counts],
                                     title='Total Counts')
                    else:
                        fig = px.bar(
                            plot_data, y=y_axis, x='count',
                            orientation='h',
                            title=f'Counts by {y_axis}',
                            labels={'count': 'Count', y_axis: y_axis}
                        )
            else:
                if x_axis == 'count' and y_axis == 'count':
                    fig.add_annotation(text="Cannot plot count vs count",
                                       xref="paper", yref="paper", x=0.5, y=0.5,
                                       showarrow=False, font_size=16)
                elif x_axis == 'count':
                    fig = px.scatter(
                        plot_data, x='count', y=y_axis,
                        color=color_dim, size=size_dim,
                        title=f'Count vs {y_axis}'
                              f"{' • color=count' if color_dim else ''}"
                              f"{' • size=count' if size_dim else ''}",
                        labels={'count': 'Count', y_axis: y_axis}
                    )
                elif y_axis == 'count':
                    if len(plot_data) < 20 and not size_dim:
                        fig = px.bar(
                            plot_data, x=x_axis, y='count',
                            title=f'Counts by {x_axis}',
                            labels={'count': 'Count', x_axis: x_axis}
                        )
                    else:
                        fig = px.scatter(
                            plot_data, x=x_axis, y='count',
                            color=color_dim, size=size_dim,
                            title=f'Counts by {x_axis}'
                                  f"{' • color=count' if color_dim else ''}"
                                  f"{' • size=count' if size_dim else ''}",
                            labels={'count': 'Count', x_axis: x_axis}
                        )
                else:
                    fig = px.scatter(
                        plot_data, x=x_axis, y=y_axis,
                        color=color_dim, size=size_dim,
                        title=f'{y_axis} vs {x_axis}'
                              f"{' • color=count' if color_dim else ''}"
                              f"{' • size=count' if size_dim else ''}",
                        labels={x_axis: x_axis, y_axis: y_axis}
                    )

        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            font=dict(size=12),
            title_font_size=16,
            showlegend=False,
            height=800
        )

        # Remove white marker borders everywhere; give a default size if not scaling by count
        fig.update_traces(
            marker=dict(
                line=dict(width=0),
            ),
            selector=dict(mode='markers')
        )

        # If not using size-by-count, ensure a readable default marker size
        if size_dim is None:
            fig.update_traces(marker=dict(size=2), selector=dict(mode='markers'))

        # Tighten margins so the plot uses more space
        fig.update_layout(
            margin=dict(l=0, r=0, t=44, b=0),  # leave a little room for the title
        )

        return fig, stats

    # Tailwind
    app.index_string = '''
    <!DOCTYPE html>
    <html>
        <head>
            {%metas%}
            <title>{%title%}</title>
            {%favicon%}
            {%css%}
            <script src="https://cdn.tailwindcss.com"></script>
        </head>
        <body>
            {%app_entry%}
            <footer>
                {%config%}
                {%scripts%}
                {%renderer%}
            </footer>
        </body>
    </html>
    '''

    app.run_server(debug=True, port=8050)