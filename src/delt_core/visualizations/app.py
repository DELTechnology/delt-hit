from collections import defaultdict

import dash
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
from dash import Dash, dcc, html, dash_table, Input, Output, State


def create_3d_figure(data: pd.DataFrame):
    if len(data.columns.drop('Count')) == 2:
        fig = px.scatter_3d(
            data,
            x='Code1', y='Code2', z='Count',
            # color='Count', size='Count',
            hover_data=['Count']
        )
    else:
        fig = px.scatter_3d(
            data,
            x='Code1', y='Code2', z='Code3',
            color='Count', size='Count',
            hover_data=['Count']
        )

    fig.update_layout(
        scene=dict(aspectmode='cube'),
        margin=dict(l=0, r=0, t=0, b=0),
        scene_camera=dict(eye=dict(x=2.0, y=2.0, z=0.75))
    )
    return fig


def parse_code_ranges(input_string, data: pd.DataFrame):
    code_ranges = defaultdict(list)
    parts = input_string.split(';')
    cols = data.columns.drop('Count')

    for col, part in zip(cols, parts):
        part = part.strip()
        subpart = part.split(',')
        for spart in subpart:
            spart = spart.strip()
            if '-' in spart:
                min_str, max_str = spart.split('-', 1)
                min_val = int(min_str.strip())
                max_val = int(max_str.strip())
                code_ranges[col].extend(range(min_val, max_val + 1))
            else:
                # single value
                if spart != '':
                    code_ranges[col].append(int(spart))

    return code_ranges


def get_default_code_ranges(data: pd.DataFrame):
    default_code_ranges = []
    for code in data.columns.drop('Count'):
        cmin = data[code].min()
        cmax = data[code].max()
        default_code_ranges.append(f"{cmin}-{cmax}")
    default_range_str = ';'.join(default_code_ranges)
    return default_range_str


def main(data: pd.DataFrame):
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    default_range_str = get_default_code_ranges(data)

    app.layout = dbc.Container(fluid=True, children=[
        # Custom header with simple styling
        html.Div(
            children=[
                html.H1(
                    "DELT Dashboard",
                    style={
                        "color": "#ffffff",
                        "margin": 0,
                        "padding": "10px 20px",
                        "fontSize": "1.5rem",
                        "fontWeight": "300",  # Thinner font weight
                        "fontFamily": "Arial, sans-serif"
                    }
                )
            ],
            style={
                "backgroundColor": "#4a90e2",  # Softer blue
                "width": "100%",
                "textAlign": "left",
                "borderRadius": "8px",  # Slight rounding for a modern feel
                "boxShadow": "0 2px 4px rgba(0, 0, 0, 0.1)",  # Subtle shadow for depth
                "marginBottom": "20px",  # Space below the header
                "marginTop": "20px"  # Space below the header
            }
        ),

        dbc.Row([
            # Left: Filters
            dbc.Col(
                dbc.Card([
                    dbc.CardHeader(html.H6("Filters", className="mb-0 fw-bold")),
                    dbc.CardBody([
                        dbc.Label("Code Ranges (e.g. '1-3;1-3;1-3'):", className="small"),
                        dcc.Input(
                            id='filter-codes',
                            type='text',
                            value=default_range_str,
                            style={'width': '100%'},
                            className="mb-3"
                        ),
                        dbc.Label("Min Count:", className="small"),
                        dcc.Input(
                            id='filter-min-count',
                            type='number',
                            value=data['Count'].min(),
                            step=1,
                            style={'width': '100%'},
                            className="mb-3"
                        ),
                        dbc.Label("Max Count:", className="small"),
                        dcc.Input(
                            id='filter-max-count',
                            type='number',
                            value=data['Count'].max(),
                            step=1,
                            style={'width': '100%'},
                            className="mb-3"
                        ),
                        dbc.Button("Filter", id='filter-button', n_clicks=0, outline=True, color="primary"),
                        dbc.Button("Reset", id='reset-button', n_clicks=0, outline=True, color="secondary", style={'marginLeft': '10px'}),
                    ])
                ], className="shadow-sm bg-white rounded"), width=3
            ),

            # Right: Data Table
            dbc.Col(
                dbc.Card([
                    dbc.CardHeader(html.H6("Data Table", className="mb-0 fw-bold")),
                    dbc.CardBody([
                        dash_table.DataTable(
                            id='data-table',
                            columns=[{"name": col, "id": col} for col in data.columns],
                            data=data.to_dict('records'),
                            page_size=10,
                            sort_action='native',
                            style_table={'overflowX': 'auto', 'width': '100%'},
                            style_cell={'textAlign': 'center', 'fontSize': '0.8rem'},
                            style_header={'fontWeight': 'bold', 'backgroundColor': '#f8f9fa'}
                        )
                    ], className="p-2")
                ], className="shadow-sm bg-white rounded"), width=9
            )
        ], className="mb-4"),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(html.H6("3D Scatter Plot", className="mb-0 fw-bold")),
                    dbc.CardBody([
                        dcc.Graph(
                            id='scatter-3d',
                            figure=create_3d_figure(data),
                            style={'width': '100%', 'height': '600px'}
                        )
                    ], className="p-2")
                ], className="shadow-sm bg-white rounded")
            ], width=12)
        ])
    ])

    @app.callback(
        Output('data-table', 'data'),
        Output('scatter-3d', 'figure'),
        Input('filter-button', 'n_clicks'),
        State('filter-codes', 'value'),
        State('filter-min-count', 'value'),
        State('filter-max-count', 'value')
    )
    def apply_filter(n_clicks, code_range_str, min_counts, max_counts):
        filtered_df = data.copy()

        if n_clicks > 0:
            code_ranges = parse_code_ranges(code_range_str, data=data)

            # Apply code range filters
            for code, vals in code_ranges.items():
                filtered_df = filtered_df[filtered_df[code].isin(vals)]

            # Apply count range filter
            if min_counts is not None:
                filtered_df = filtered_df[filtered_df['Count'] >= min_counts]
            if max_counts is not None:
                filtered_df = filtered_df[filtered_df['Count'] <= max_counts]

        table_data = filtered_df.to_dict('records')
        fig = create_3d_figure(filtered_df)
        return table_data, fig

    @app.callback(
        Output('filter-codes', 'value'),
        Output('filter-min-count', 'value'),
        Output('filter-max-count', 'value'),
        Input('reset-button', 'n_clicks')
    )
    def reset_filters(n_clicks):
        # Reset to default values when "Reset" is clicked
        if n_clicks > 0:
            cmin = data['Count'].min()
            cmax = data['Count'].max()
            return default_range_str, cmin, cmax
        return dash.no_update, dash.no_update, dash.no_update

    app.run_server(debug=True)


if __name__ == '__main__':
    # Sample DataFrame
    data = pd.DataFrame({
        'Code1': [1, 2, 1, 3, 2, 3, 1],
        'Code2': [1, 1, 2, 3, 2, 1, 3],
        'Code3': [1, 2, 3, 1, 3, 2, 2],
        'Count': [230, 150, 400, 500, 100, 350, 210]
    })

    main(data=data)
