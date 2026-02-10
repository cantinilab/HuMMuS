FROM ubuntu:24.04

RUN apt update -y

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata
RUN apt install -y build-essential  \
                   gcc \
                   graphviz \
                   cmake \
                   git \
                   wget \
                   curl \
                   vim \
                   bzip2 \
                   libbz2-dev \
                   file \
                   libxrender-dev \
                   libxext-dev \
                   postgresql \
                   postgresql-contrib \
                   swig \
                   python3-dev
RUN apt -y install libudunits2-dev libgdal-dev libgeos-dev libproj-dev 
RUN apt -y install zlib1g-dev libabsl-dev

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
RUN cp /opt/HuMMuS/hummuspy/src/hummuspy/utils/logging/logging_conf.ini /opt/miniconda3/envs/rhummus_env/lib/python3.10/site-packages/hummuspy/utils/logging/
RUN Rscript install_rdeps.R

ENV PATH=/opt/miniconda3/bin:$PATH
