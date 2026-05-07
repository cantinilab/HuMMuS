# Perform analysis on the multilayer structure
rule analysis_hummus:
    input:
        marker_step1 = os.path.join( config['working_dir'], "ready", "analysis_step2_done.txt")
    output:
        marker_step2 = os.path.join( config['working_dir'], "ready", "analysis_done.txt")
    params:
        working_dir = config["working_dir"],
        config_file = config["configuration_file"]
    shell:
        "python scripts/exec_analysis_hummus.py {params.working_dir} {params.config_file}"