FROM conda/miniconda3-centos7

LABEL gitUrl="https://gitlab.curie.fr/data-analysis/atac-seq.git"
LABEL gitCommit="083a87393264c306352457d9271639f4fa43da06 / devel"

# real path from baseDir: /home/nservant/Apps/geniac/build/workDir/recipes/conda/bowtie2.yml
ADD bowtie2.yml /opt/bowtie2.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/bowtie2.yml \
&& echo "source activate bowtie2-2.3.5" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/bowtie2-2.3.5/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
