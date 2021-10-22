FROM conda/miniconda3-centos7

LABEL gitUrl="https://gitlab.curie.fr/data-analysis/atac-seq.git"
LABEL gitCommit="083a87393264c306352457d9271639f4fa43da06 / devel"

# real path from baseDir: /home/nservant/Apps/geniac/build/workDir/recipes/conda/bedtobigbed.yml
ADD bedtobigbed.yml /opt/bedtobigbed.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/bedtobigbed.yml \
&& echo "source activate bedtobigbed-377" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/bedtobigbed-377/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
