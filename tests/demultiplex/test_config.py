from delt_core.demultiplex.utils import Config


def test_save_load_config(tmp_path):
    config_data = {
        "Root": "/path/to/root",
        "Experiment": {"Name": ''},
        "Selection": {
            "SelectionFile": "/path/to/selection_file",
            "FASTQFile": "/path/to/fastq_file",
            "Library": "/path/to/library",
        },
        "Structure": {
            "S1": {
                "MaxErrorRate": 0.0,
                "Indels": False
            }
        }
    }
    config = Config(**config_data)
    yaml_path = tmp_path / "config.yml"
    yaml_path = config.write_yaml(yaml_path)
    loaded_config = Config.from_yaml(yaml_path).to_dict()
    loaded_config.pop('Simulation')
    assert config_data == loaded_config
