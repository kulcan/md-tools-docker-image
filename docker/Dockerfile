FROM ubuntu:22.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt upgrade -y && \
    apt-get install -y build-essential wget bzip2 ca-certificates curl git nano && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget http://vergil.chemistry.gatech.edu/psicode-download/Psi4conda-1.5-py38-Linux-x86_64.sh -O ~/miniconda.sh

RUN /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN pip install ase

RUN conda install -c conda-forge mamba
RUN apt-get update
RUN apt-get install -y libgmp3-dev
RUN mamba install -y sage -c conda-forge

RUN conda install -c conda-forge jupyterlab
RUN conda install -c conda-forge openmm

COPY 1ubq.pdb ./
COPY test.ipynb ./

EXPOSE 8888
ENTRYPOINT ["jupyter", "notebook", "--no-browser", "--ip=0.0.0.0","--NotebookApp.token=''", "--NotebookApp.password=''", "--allow-root"]
