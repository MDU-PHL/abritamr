FROM continuumio/miniconda3:4.11.0 as build

ARG amrfinder_version=3.10.36

RUN apt-get --yes update && apt-get --yes upgrade && apt-get install --yes curl wget gzip unzip git && apt-get clean

RUN conda update conda

RUN conda install -c conda-forge -c bioconda -c defaults file blast=2.9=pl526h3066fca_4 ncbi-amrfinderplus=$amrfinder_version

RUN amrfinder -u

RUN cd /opt && mkdir -p abritamr

COPY . /opt/abritamr 

RUN cd /opt/abritamr && pip install .

RUN conda remove --force curl entrez-direct tqdm idna pip setuptools wheel asn1crypto ca-certificates cryptography expat chardet libcurl libedit libssh2 ncurses openssl conda-package-handling readline pysocks pycosats pycparser pyopenssl sqlite tk urllib3 requests yaml ruamel_yaml

RUN cd /opt/conda && rm -rf pkgs && rm -rf man/man*

FROM bitnami/minideb:buster

COPY --from=build /opt/conda /opt/conda

ENV PATH=/opt/conda/bin:$PATH

CMD [ "/bin/bash" ]
