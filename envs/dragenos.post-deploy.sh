#!/usr/bin/env bash
set -e



mkdir -p ${CONDA_PREFIX}/software
cd ${CONDA_PREFIX}/software

git clone https://github.com/populationgenomics/DRAGMAP.git

cd ${CONDA_PREFIX}/software/DRAGMAP

HAS_GTEST=0 make

cp ${CONDA_PREFIX}/software/DRAGMAP/build/release/dragen-os ${CONDA_PREFIX}/bin


