from scripts.pca import *
import os
import pandas as pd
from sklearn.preprocessing import scale


def samplefile_stat_path(samplefile, suffix):
    return pj(get_samplefile_folder(samplefile), f"{samplefile}.{suffix}")

PCA_DONE_FILES = expand('PCA_done_{samplefile}', samplefile=SAMPLE_FILES)
PCA_HDF5_FILES = [samplefile_stat_path(samplefile, "coverage.hdf5") for samplefile in SAMPLE_FILES]


rule PCA_all:
    input: PCA_DONE_FILES
    default_target: True



rule clustering_samples:
    input: PCA_HDF5_FILES
    output: stat = pj(STAT, 'PCA_stat.txt'),
            file_done = PCA_DONE_FILES
            # groups = pj(STAT, 'PCA_groups.tsv')
    conda: CONDA_PCA
    resources:
        n=8,
        mem_mb=32000
    run:
        os.makedirs(os.path.dirname(str(output.stat)), exist_ok=True)
        hdf5_files = list(input)
        PCA_t, expl_var = perform_ipca(hdf5_files=hdf5_files)
        # best_dict = get_scores_and_labels(combinations,data=PCA_t)
        PCA_2 = PCA_t[[0, 1]]
        data_2 = scale(PCA_2)
        PCA_3 = PCA_t[[0, 1, 2]]
        data_3 = scale(PCA_3)
        PCA_4 = PCA_t[[0, 1, 2, 3]]
        data_4 = scale(PCA_4)
        PCA_6 = PCA_t[[0, 1, 2, 3, 4, 5]]
        data_6 = scale(PCA_6)
        PCA_10 = PCA_t
        data_10 = scale(PCA_10)
        best_dict_2 = get_scores_and_labels(combinations,data=data_2)
        best_dict_3 = get_scores_and_labels(combinations,data=data_3)
        best_dict_4 = get_scores_and_labels(combinations,data=data_4)
        best_dict_6 = get_scores_and_labels(combinations,data=data_6)
        best_dict_10 = get_scores_and_labels(combinations,data=data_10)
        PCA_t['cluster_2_cor_cov'] = best_dict_2['best_labels']
        PCA_t['cluster_3_cor_cov'] = best_dict_3['best_labels']
        PCA_t['cluster_4_cor_cov'] = best_dict_4['best_labels']
        PCA_t['cluster_6_cor_cov'] = best_dict_6['best_labels']
        PCA_t['cluster_10_cor_cov'] = best_dict_10['best_labels']
        clusters = ['cluster_2_cor_cov', 'cluster_3_cor_cov', 'cluster_4_cor_cov', 'cluster_6_cor_cov',
                    'cluster_10_cor_cov']
        with open(output.stat, 'w') as f:
            print(f"Cumulative variance = {expl_var}",file=f)
            print(best_dict_2, file=f)
            print(best_dict_3,file=f)
            print(best_dict_4,file=f)
            print(best_dict_6,file=f)
            print(best_dict_10,file=f)
        for cluster in clusters:
            for hdf5_file in hdf5_files:
                add_cluster_info_to_hdf5(hdf5_file, clustred_data=PCA_t, cluster=cluster)
        for done_file in output.file_done:
            with open(done_file, 'w') as f:
                f.write("")
