# Package `delt-core`
Core functionalities to work with DECL libraries

# Installation

Install the package for development purposes.
Navigate to the root folder of the package and run the following command:

```bash
pip install -e ".[dev,test]"
```

Create a `.env` file to store environment configurations, keys and secrets.
```bash
touch .env
```
These configurations can be accessed using the `python-dotenv` package.


# CLI

In a terminal, run the following commands to test the CLI:

```bash
delt-core --help
delt-core welcome # Welcome to DEL Technology!
delt-core process # Processing...
```