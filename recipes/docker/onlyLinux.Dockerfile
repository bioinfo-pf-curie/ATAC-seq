FROM centos:7

LABEL gitUrl="https://gitlab.curie.fr/data-analysis/atac-seq.git"
LABEL gitCommit="083a87393264c306352457d9271639f4fa43da06 / devel"

RUN yum install -y which \
&& yum clean all

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
