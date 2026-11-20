#!/usr/bin/env bash
set -e

conda install -y conda-forge::unzip conda-forge::wget

mkdir -p ${CONDA_PREFIX}/share
cd ${CONDA_PREFIX}/share/

wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
rm ${CONDA_PREFIX}/bin/gatk
ln -s ${CONDA_PREFIX}/share/gatk-4.5.0.0/gatk ${CONDA_PREFIX}/bin/gatk

conda remove -y unzip wget
