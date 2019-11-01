Bootstrap: docker
From: continuumio/miniconda3:4.5.12

%help
A Singularity image for AMR finder plus 

%labels
Maintainer Kristy Horan
Build 1.0
amrfinderplus-3.1.1b
db 2019-09-06.1

%environment
export PATH=/opt/conda/bin:$PATH
%post
 # set versions of software to install
  export PATH=/opt/conda/bin:$PATH 
  conda --version
  conda config --add channels conda-forge
  conda config --add channels defaults
  conda config --add channels r
  conda config --add channels bioconda
 
# install amrfinderplus

  conda install -y -c bioconda ncbi-amrfinderplus
  conda update -c bioconda -y ncbi-amrfinderplus
  

# get test files
  curl -O https://raw.githubusercontent.com/ncbi/amr/v3b/test_dna.fa \
  -O https://raw.githubusercontent.com/ncbi/amr/v3b/test_prot.fa \
  -O https://raw.githubusercontent.com/ncbi/amr/v3b/test_prot.gff \
  -O https://raw.githubusercontent.com/ncbi/amr/v3b/test_both.expected \
  -O https://raw.githubusercontent.com/ncbi/amr/v3b/test_dna.expected \
  -O https://raw.githubusercontent.com/ncbi/amr/v3b/test_prot.expected 
 
  
  amrfinder -u
  echo "Done"

%test
	amrfinder --plus -n test_dna.fa -O Campylobacter > test_dna.got
	diff test_dna.expected test_dna.got