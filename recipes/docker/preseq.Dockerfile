FROM conda/miniconda3-centos7

LABEL gitUrl="https://gitlab.curie.fr/data-analysis/atac-seq.git"
LABEL gitCommit="083a87393264c306352457d9271639f4fa43da06 / devel"

# real path from baseDir: /home/nservant/Apps/geniac/build/workDir/recipes/conda/preseq.yml
ADD preseq.yml /opt/preseq.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/preseq.yml \
&& echo "source activate preseq-2.0.3" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/preseq-2.0.3/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
