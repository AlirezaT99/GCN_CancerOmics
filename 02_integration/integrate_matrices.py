__author__ = "Alireza Tajmirriahi"
__version__ = "1.0.0"
__maintainer__ = "Alireza Tajmirriahi"
__email__ = "alireza.tajmirriahi@sharif.edu"
__project__ = "CancerGenomics"

import argparse
import os

import numpy as np
import pandas as pd

from tqdm import tqdm

DATA_DIR = '/PROJECTS/Taj/1_PreprocessData/data'
MATRICES_PATH = '/PROJECTS/Taj/1_PreprocessData/output'
OUT_DIR = 'ICGC'
MIN_SAMPLES = 10  # 50
TRAIN_SPLIT = 0.70


def aggregate_dataframes(data):
    cols = sorted(set.union(*[set(df.columns) for df in data.values()]))
    to_merge = []
    for label, df in data.items():
        new_df = df.reindex(columns=cols, fill_value=0)
        new_df['label'] = label
        to_merge.append(new_df)
    return pd.concat(to_merge)


def agg_mutation(mat, cancer_type, file_name):
    df = pd.read_csv(f'{args.DATA_DIR}/{cancer_type}/{file_name}', sep='\t')
    mutation_count = df.groupby(['icgc_donor_id', 'gene_symbol']).size()

    result = mat.copy()
    for row in tqdm(result.index):
        cols = mat.columns[np.where(mat.loc[row] == 1)]
        for col in cols:
            result.loc[row][col] = mutation_count[row][col]
    return result


def agg_expression(mat, cancer_type, file_name):
    df = pd.read_csv(f'{args.DATA_DIR}/{cancer_type}/{file_name}')
    df = df.drop(columns=['Unnamed: 0'], axis=1)
    df = df.set_index(['icgc_donor_id', 'gene_id'])
    df = df[~df.index.duplicated(keep='first')]

    result = mat.copy().astype(np.float64)
    for row in tqdm(result.index):
        cols = mat.columns[np.where(mat.loc[row] == 1)]
        for col in cols:
            result.loc[row][col] = df.loc[row, col].normalized_expression_value
    return result


def main():
    file_paths = [*os.walk(MATRICES_PATH)][0][2]
    all_dfs = list(map(lambda p: pd.read_csv(f'{args.MATRICES_PATH}/{p}', index_col=0, delimiter='\t'), file_paths))
    all_labels = list(map(lambda p: p[7:-4], file_paths))

    omic_paths = ['symbol_mutation.tsv', 'expression_data.tsv']
    mut_data, exp_data = dict(), dict()

    for label, df in zip(all_labels, all_dfs):
        if len(df) >= args.MIN_SAMPLES:
            print('Aggregating', label, '...')
            mut_data[label] = agg_mutation(df, label, omic_paths[0])
            exp_data[label] = agg_expression(df, label, omic_paths[1])

    print('Merging dataframes', end=' ')
    merged_mut = aggregate_dataframes(mut_data)
    merged_exp = aggregate_dataframes(exp_data)
    print(f'Done')

    assert merged_mut.shape == merged_exp.shape
    labels = merged_mut.label

    train_indices = np.random.choice([True, False], len(merged_mut), p=[args.TRAIN_SPLIT, 1 - args.TRAIN_SPLIT])

    labels_train = labels[train_indices]
    labels_test = labels[~train_indices]
    mut_train = merged_mut[train_indices]
    mut_test = merged_mut[~train_indices]
    exp_train = merged_exp[train_indices]
    exp_test = merged_exp[~train_indices]

    if not os.path.exists(args.OUT_DIR):
        os.mkdir(args.OUT_DIR)
    labels_train.to_csv(f'{args.OUT_DIR}/labels_tr.csv', index=False, header=False)
    labels_test.to_csv(f'{args.OUT_DIR}/labels_te.csv', index=False, header=False)
    merged_mut.columns.to_frame().to_csv(f'{args.OUT_DIR}/1_featname.csv', index=False, header=False)
    merged_exp.columns.to_frame().to_csv(f'{args.OUT_DIR}/2_featname.csv', index=False, header=False)
    mut_train.to_csv(f'{args.OUT_DIR}/1_tr.csv', index=False, header=False)
    mut_test.to_csv(f'{args.OUT_DIR}/1_te.csv', index=False, header=False)
    exp_train.to_csv(f'{args.OUT_DIR}/2_tr.csv', index=False, header=False)
    exp_test.to_csv(f'{args.OUT_DIR}/2_te.csv', index=False, header=False)

    print('num classes=', len(mut_data.keys()))
    # run_mogonet(num_class=len(mut_data.keys()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--data-dir', action='store', dest='DATA_DIR', type=str,
                        default='/PROJECTS/Taj/1_PreprocessData/data', help='Path to data')
    parser.add_argument('--matrices-path', action='store', dest='MATRICES_PATH', type=str,
                        default='/PROJECTS/Taj/1_PreprocessData/output', help='Path to preprocessed matrices')
    parser.add_argument('--out-dir', action='store', dest='OUT_DIR', type=str, default='ICGC',
                        help='The output directory')
    parser.add_argument('--min-samples', action='store', dest='MIN_SAMPLES', type=int, default=10,
                        help='Minimum samples required to include a cancer type')
    parser.add_argument('--train-split', action='store', dest='TRAIN_SPLIT', type=float, default=0.70,
                        help='train/test split. By default set to 0.70')

    args = parser.parse_args()
    main()
