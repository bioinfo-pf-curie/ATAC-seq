FROM conda/miniconda3-centos7

LABEL gitUrl="https://gitlab.curie.fr/data-analysis/atac-seq.git"
LABEL gitCommit="083a87393264c306352457d9271639f4fa43da06 / devel"

# real path from baseDir: /home/nservant/Apps/geniac/build/workDir/recipes/conda/deeptools.yml
ADD deeptools.yml /opt/deeptools.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/deeptools.yml \
&& echo "source activate deeptools-3.2.1" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/deeptools-3.2.1/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
