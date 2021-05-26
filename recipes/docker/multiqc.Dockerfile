FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/atac-seq.git"
LABEL gitCommit="be7438449678fae2a43874ac6ec57616eaa801f9 / devel"

# real path from baseDir: /home/nservant/Apps/geniac/build/workDir/recipes/conda/multiqc.yml
ADD multiqc.yml /opt/multiqc.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/multiqc.yml \
&& echo "source activate multiqc-1.9" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/multiqc-1.9/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
