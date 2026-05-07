# generate layers and networks
rule generate_networks_layers:
    input:
        marker_step1 = os.path.join( config['working_dir'], "ready", "analysis_step1_done.txt")
    output:
        marker_step2 = os.path.join( config['working_dir'], "ready", "analysis_step2_done.txt")
    params:
        working_dir = config["working_dir"],
        config_file = config["configuration_file"]
    shell:
        "python scripts/exec_generate_networks_layers.py {params.working_dir} {params.config_file}"