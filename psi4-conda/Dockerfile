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

ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

COPY sample-psi4.py ./

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
