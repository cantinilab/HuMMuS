# preprocess assays data
rule process_data:
    output:
        marker_step1 = os.path.join( config['working_dir'], "ready", "analysis_step1_done.txt")
    params:
        working_dir = config["working_dir"],
        config_file = config["configuration_file"]
    shell:
        "python scripts/exec_process_data.py {params.working_dir} {params.config_file}"