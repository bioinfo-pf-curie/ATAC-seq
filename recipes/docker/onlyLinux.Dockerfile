FROM centos:7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/atac-seq.git"
LABEL gitCommit="be7438449678fae2a43874ac6ec57616eaa801f9 / devel"

RUN yum install -y which \
&& yum clean all

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
