FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/atac-seq.git"
LABEL gitCommit="be7438449678fae2a43874ac6ec57616eaa801f9 / devel"

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
