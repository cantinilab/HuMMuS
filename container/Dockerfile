FROM ubuntu:20.04

RUN mkdir /opt
WORKDIR /opt
COPY ./environment.yml .
COPY ./install_rdeps.R .

RUN mkdir -p /opt/miniconda3
RUN cd /opt/miniconda3 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /opt/miniconda3 -b -f
RUN rm Miniconda3-latest-Linux-x86_64.sh
RUN export PATH=/opt/miniconda3/bin:$PATH
RUN export conda_NO_LOW_SPEED_LIMIT=1   # no timeout on conda commands
RUN export PIPENV_INSTALL_TIMEOUT=3600   # 1h timeout on pip install commands 
RUN export PIP_DEFAULT_TIMEOUT=3600

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda update conda -y
RUN . /opt/miniconda3/etc/profile.d/conda.sh

RUN conda config --add channels conda-forge
RUN conda env create --name rhummus_env --file=/opt/environment.yml

RUN conda activate rhummus_env
RUN cd /opt
RUN git clone https://github.com/YasCoMa/HuMMuS
    
RUN cd /opt/HuMMuS
RUN git switch task_improve_parameter_doc

RUN cd /opt
RUN pip install /opt/HuMMuS/hummuspy
RUN Rscript install_rdeps.R

ENV PATH=/opt/miniconda3/bin:$PATH
