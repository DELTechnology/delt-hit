# AGENTS.md — project-level guidance for AI coding agents

## Python environment

This project uses **pixi** for environment and dependency management. Always run Python via:

```bash
pixi run python <script.py>
pixi run python -c "<snippet>"
pixi run pytest ...
```

Never use `uv`, `uv run`, `pip install`, or bare `python3` to execute project code. If you need to run a one-liner, use `pixi run python -c "..."`.
