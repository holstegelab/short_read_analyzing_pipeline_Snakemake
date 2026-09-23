## Building the Apptainer environment

**Run all commands from the directory where this file is located**

The workflow uses a patched Snakemake branch with support for containerizing post deploy scripts in Conda environments:

    pip install git+https://github.com/maarten-k/snakemake@containerize_postdeploy

Generate the container definition from the dedicated environment-generation workflow. No aditional settings are needed to run this:

    snakemake --conda-create-envs-only --containerize apptainer -s envs/make_all_envs.smk > short_read_env.def


Build the Apptainer image:

    apptainer build -B $PWD:/tmp/bind short_read_env.sif short_read_env.def

### `gatk_cnv` environment

The `gatk_cnv` environment is currently not actively used. Its post-processing requirements appear to be replaceable by installing `gcnvkernel` directly:

    conda install bioconda::gcnvkernel

This should be verified against the existing `gatk_cnv` post-processing workflow before removing the environment entirely.