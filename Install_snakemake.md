1. Download and install zslurm according to instructions: https://github.com/holstegelab/zslurm  (NB! install all dependencies, without them there will be no error, but whole thing won't work)
2. Download and install https://github.com/holstegelab/snakemake-executor-plugin-zslurm (poetry install)
3. Download snakemake https://github.com/mhulsman/snakemake 
   switch to development_v3 branch
   python -m pip install -e .

Caveats:
* previously installed snakemake could cause troubles even if deleted from conda
* installing or updating plugin can cause updating snakemake to the wrong version
* 