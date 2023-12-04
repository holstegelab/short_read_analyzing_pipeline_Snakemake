from scripts.pca import *
import pandas as pd
from sklearn.preprocessing import scale


rule PCA_all:
    input: expand('PCA_done_{samplefile}', samplefile = SAMPLE_FILES)
    default_target: True



rule clustering_samples:
    input: expand('{samplefile}.coverage.hdf5',  samplefile = SAMPLE_FILES)
    output: stat = pj(STAT, 'PCA_stat.txt'),
            file_done = touch(temp('PCA_done_{samplefile}'))
            # groups = pj(STAT, 'PCA_groups.tsv')
    conda: CONDA_PCA
    run:
        shell("touch {output.stat}")
        hdf5_files = input
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
        with open(pj(STAT,'PCA_stat.txt'),'a') as f:
            print(f"Cumulative variance = {expl_var}",file=f)
            print(best_dict_2, file=f)
            print(best_dict_3,file=f)
            print(best_dict_4,file=f)
            print(best_dict_6,file=f)
            print(best_dict_10,file=f)
        for cluster in clusters:
            for SF in SAMPLE_FILES:
                add_cluster_info_to_hdf5(f'{SF}', clustred_data=PCA_t, cluster=cluster)