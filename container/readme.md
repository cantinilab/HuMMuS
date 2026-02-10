## Usage of container for HuMMuS

Basically, you have two choices of virtual machine engines: singularity or docker

#### Singularity
- Build image:
    `bash setup_vm_singularity.sh`
- Run image entering in interactive shell:
    `singularity shell ~/hummus_vm/hummus.simg`
- Run image executing an analysis script:
    `singularity exec --nv --cleanenv ~/hummus_vm/hummus.simg python analysis.py`
    
#### Docker
- Build image:
    `docker build -t cantinilab/hummusvm .`
- Run image entering in interactive shell:
    `sudo docker run -i -t cantinilab/hummusvm --entrypoint /bin/bash`
- Run image executing an analysis script:
    `sudo docker exec -d cantinilab/hummusvm python analysis.py`

