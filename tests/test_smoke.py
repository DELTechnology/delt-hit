"""Smoke test to verify the package installs and CLI loads."""
import subprocess
import sys


def test_import():
    """Test that the main package can be imported."""
    import delt_hit

    assert delt_hit is not None


def test_cli_help():
    """Test that the CLI --help flag works."""
    result = subprocess.run(
        [sys.executable, "-m", "delt_hit.cli.main", "--help"],
        capture_output=True,
        text=True,
    )
    # jsonargparse exits 0 on --help and prints to stdout
    assert result.returncode == 0
    output = result.stdout + result.stderr
    assert "usage" in output.lower()
