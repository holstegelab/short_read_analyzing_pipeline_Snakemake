import numpy as np
import os
import pandas as pd
from sklearn.decomposition import IncrementalPCA
from sklearn.cluster import DBSCAN
# from sklearn.model_selection import ParameterGrid
from sklearn.preprocessing import scale
from sklearn.metrics import silhouette_score as ss
import pickle
import h5py
from sklearn.preprocessing import scale
import read_samples
from common import *
import utils

group_name = 'coverage'
#
main_chrs = ['chr1:', 'chr2:', 'chr3:',
             'chr4:', 'chr5:', 'chr6:',
             'chr7:', 'chr8:', 'chr9:',
             'chr10:', 'chr11:', 'chr12:',
             'chr13:', 'chr14:', 'chr15:',
             'chr16:', 'chr17:', 'chr18:',
             'chr19:', 'chr20:', 'chr21:',
             'chr22:', 'chrX:', 'chrY:']
#
#
# def load_hdf5_data(sample_file, df, group_name):
#     with h5py.File(sample_file, 'r') as f:
#         data_group = f[group_name]
#         samples = [s.decode('utf-8') for s in data_group['samples'][()]]
#         chrs = data_group['chrom'][()]
#         start = data_group['start'][()]
#         end = data_group['end'][()]
#         cov = data_group['coverage'][()]
#         index = [f"{chrom.decode('utf-8')}:{s}-{e}" for chrom, s, e in zip(chrs, start, end)]
#         # Create a temporary DataFrame for the new data
#         temp_df = pd.DataFrame(
#             columns=samples,
#             data=cov.T,
#             index=index
#         )
#         df_filtered = temp_df[temp_df.index.str.startswith(tuple(main_chrs))]
#         df = pd.concat([df, df_filtered], axis=1)
#
#     return df, len(samples)
#
# def perform_pca(df, n_comp):
#     pca = PCA(n_components = n_comp)
#     pca.fit(df.T)
#     transformed_data = pca.transform(df.T)
#     scores = pd.DataFrame(transformed_data, index=df.T.index)
#     explained_variance_ratio = pca.explained_variance_ratio_
#     cum_expl_var = np.cumsum(explained_variance_ratio)
#     return scores, cum_expl_var
# MAX_DELTA = np.inf
# previous_components = 0
def perform_ipca(hdf5_files, batch_size = 1000, max_delta = 1e-5, ncomponents = 10):
    total_intervals = 0
    for file in hdf5_files:
        with h5py.File(file, 'r') as f:
            total_intervals = max(total_intervals, f[group_name]['coverage'].shape[1])

    ipca = IncrementalPCA(n_components=ncomponents)
    MAX_DELTA = np.inf
    previous_components = None
    i = 0
    sample_names = []
    while i < total_intervals and MAX_DELTA > max_delta:
        df = pd.DataFrame()
        batch_sample_names = []
        for file in hdf5_files:
            with h5py.File(file, 'r') as f:
                data_group = f[group_name]
                samples = [s.decode('utf-8') for s in data_group['samples'][()]]
                chrs = data_group['chrom'][i:i + batch_size]
                start = data_group['start'][i:i + batch_size]
                end = data_group['end'][i:i + batch_size]
                index = [f"{chrom.decode('utf-8')}:{s}-{e}" for chrom, s, e in zip(chrs, start, end)]
                bases_mapped_raw = np.asarray(data_group['bases_mapped'][:], dtype=np.float64)
                valid_samples = np.isfinite(bases_mapped_raw) & (bases_mapped_raw > 0)
                if not np.all(valid_samples):
                    skipped = [sample for sample, valid in zip(samples, valid_samples) if not valid]
                    print(f"Skipping samples with invalid bases_mapped in {file}: {', '.join(skipped)}")
                if not np.any(valid_samples):
                    continue
                samples = [sample for sample, valid in zip(samples, valid_samples) if valid]
                batch_sample_names.extend(samples)
                cov = np.asarray(data_group['coverage'][valid_samples, i:i + batch_size], dtype=np.float64)
                bases_mapped = bases_mapped_raw[valid_samples][:, np.newaxis] / 1000000000.0
                corrected_cov = cov / bases_mapped
                temp_df = pd.DataFrame(
                            columns=samples,
                            data=corrected_cov.T,
                            index=index
                        )
                temp_filtred = temp_df[temp_df.index.str.startswith(tuple(main_chrs))]
                df = pd.concat([df, temp_filtred], axis=1)

        df = df.replace([np.inf, -np.inf], np.nan)
        nonfinite_rows = df.isna().any(axis=1)
        if nonfinite_rows.any():
            print(f"Skipping {int(nonfinite_rows.sum())} intervals with non-finite normalized coverage in PCA batch starting at interval {i}")
            df = df.loc[~nonfinite_rows]
        if df.empty:
            i = i + batch_size
            continue
        if len(sample_names) == 0:
            sample_names = batch_sample_names
        ipca.partial_fit(df)
        components = ipca.components_
        explained_variance_ratio = ipca.explained_variance_ratio_
        cum_expl_var = np.cumsum(explained_variance_ratio)
        print(f"Cumulative variance = {cum_expl_var}")
        if previous_components is not None:
            delta = np.abs(components - previous_components)
            MAX_DELTA = np.max(delta)
        i = i + batch_size
        previous_components = components
    if previous_components is None:
        raise ValueError("No finite normalized coverage rows were available for PCA.")
    pca_components = pd.DataFrame(ipca.components_.T, index=sample_names)
    explained_variance_ratio = ipca.explained_variance_ratio_
    cum_expl_var = np.cumsum(explained_variance_ratio)
    return pca_components, cum_expl_var


epsilons = np.linspace(0.01, 1, num = 30)

min_samples = np.arange(2, 39, step = 3)
combinations = list(itertools.product(epsilons, min_samples))
N = len(combinations)


def get_scores_and_labels(combinations, data):
    scores = []
    all_labels_list = []
    X = scale(data)
    fallback_labels = np.zeros(X.shape[0], dtype=int)
    for i, (eps, num_samples) in enumerate(combinations):
        dbscan_cluster_model = DBSCAN(eps=eps, min_samples=num_samples).fit(X)
        labels = dbscan_cluster_model.labels_
        labels_set = set(labels)
        num_clusters = len(labels_set)
        if -1 in labels_set:
            num_clusters -= 1
        if (num_clusters < 2):
            scores.append(-10)
            all_labels_list.append(fallback_labels)
            c = (eps, num_samples)
            print(f"Combination {c} on iteration {i + 1} of {N} has {num_clusters} clusters. Moving on")
            continue
        scores.append(ss(X, labels))
        all_labels_list.append(labels)
        print(f"Index: {i}, Score: {scores[-1]}, NumClusters: {num_clusters}")

    best_index = np.argmax(scores)
    best_parameters = combinations[best_index]
    best_labels = all_labels_list[best_index]
    best_score = scores[best_index]
    if best_score == -10:
        return {'best_epsilon': None,
                'best_min_samples': None,
                'best_labels': fallback_labels,
                'best_score': None,
                'fallback': 'single_cluster'}
    return {'best_epsilon': best_parameters[0],
            'best_min_samples': best_parameters[1],
            'best_labels': best_labels,
            'best_score': best_score,
            'fallback': None}

def add_cluster_info_to_hdf5(SF, clustred_data, cluster):
    file_path = SF if os.path.exists(SF) else f'{SF}.coverage.hdf5'
    with h5py.File(file_path, 'a') as f:
        data_group = f['coverage']

        # Extract sample names from the 'coverage' group
        samples_in_file = [s.decode('utf-8') for s in data_group['samples'][()]]
        # print(samples_in_file)
        mask = [sample in samples_in_file for sample in clustred_data.index]
        # print(mask)
        # Missing samples were excluded from PCA and get DBSCAN noise label -1.
        cluster_values = pd.Series(clustred_data[cluster].values, index=clustred_data.index)
        relevant_data = np.array([cluster_values.get(sample, -1) for sample in samples_in_file], dtype=int)
        # print(len(relevant_data))
        # print(len(list(data_group['samples'])))

        if cluster in data_group:
            # If 'group' dataset exists, modify its values with the relevant data
            data_group[cluster][:] = relevant_data
        else:
            # If 'group' dataset doesn't exist, create a new one and write the relevant data into it
            data_group.create_dataset(cluster, data=relevant_data)
        # Print samples and corresponding 'group' values
        samples = [s.decode('utf-8') for s in data_group['samples'][()]]
        group_values = data_group[cluster][()]

        for sample, group_value in zip(samples, group_values):
            return f"Sample: {sample}, Group: {group_value}"

# best_dict = get_scores_and_labels(combinations, X = data)
#
# print(best_dict)
# PCA_t['cluster'] = best_dict['best_labels']
