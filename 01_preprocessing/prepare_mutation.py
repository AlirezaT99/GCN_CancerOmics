import os
import argparse

import numpy as np
import pandas as pd

from tqdm import tqdm

tqdm.pandas()

DATA_DIR = './data'


def read_data(folder: str, file_path: str):
    if file_path.endswith('.tsv'):
        return pd.read_csv(f'{DATA_DIR}/{folder}/{file_path}', sep='\t', header=0, low_memory=False)
    return None


def read_data_csv(folder: str, file_path: str):
    return pd.read_csv(f'{DATA_DIR}/{folder}/{file_path}', header=0, low_memory=False)


def read_chunk_by_chunk(folder: str, file_path: str, columns=None):
    df = pd.DataFrame()
    for chunk in pd.read_csv(f"{DATA_DIR}/{folder}/{file_path}", sep='\t', header=0, low_memory=False, chunksize=1e6):
        df = pd.concat([df, chunk[columns] if columns else chunk], ignore_index=True)
    return df


def save_tsv(df, output_path, file_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    df.to_csv(f'{output_path}/{file_path}', sep='\t')


def get_mutation_data(cancer_type, mutation_path, mutation_type=None):
    columns = ['icgc_donor_id', 'gene_affected', 'mutation_type']
    data = read_chunk_by_chunk(cancer_type, mutation_path, columns)
    if mutation_type:
        data = data[data['mutation_type'] == mutation_type] \
            .drop(columns=['mutation_type'])
    return data.dropna()


def get_genes(genes_path, gene_class=None):
    genes = pd.read_csv(genes_path, sep='\t', header=0)
    genes.gene_symbol = list(map(lambda g: g[1:-1], genes.gene_symbol))
    if gene_class:
        genes = genes[genes['gene_class'] == gene_class]
    genes = genes[['gene_name', 'gene_symbol']] \
        .rename({'gene_name': 'gene_ensembl_id'}, axis=1)
    return genes


def perform_analysis(args, cancer_type):
    genes = get_genes(args.genes_path)
    print('Converting', cancer_type, end='...')


    ### Mutation
    mut = get_mutation_data(cancer_type, args.mutation_path, mutation_type='single base substitution')
    mut_data = mut.rename({'gene_affected': 'gene_ensembl_id'}, axis=1)
    sign_mut_samples = pd.merge(genes, mut_data, how='left', on='gene_ensembl_id') \
        .drop(columns=['gene_ensembl_id']) \
        .drop_duplicates() \
        .dropna()
    sign_mut_samples.to_csv(f'{DATA_DIR}/{cancer_type}/symbol_mutation.tsv', sep='\t')
    print('done')


def run(args):
    if not args.cancer_type:
        if args.run_all:
            sub_folders = [f.name for f in os.scandir(args.data_path) if f.is_dir()]
            for cancer_type in sub_folders:
                perform_analysis(args, cancer_type)
        else:
            raise Exception('Either set --cancer-type or set run_all to True')
    if not os.path.exists(f'{args.data_path}/{args.cancer_type}'):
        raise Exception('arg --cancer-type is not a valid directory')
    perform_analysis(args, args.cancer_type)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--cancer-type', type=str)  # , default='Test'
    parser.add_argument('--run-all', type=bool, default=False)

    parser.add_argument('--data-path', type=str, default='./data')
    parser.add_argument('--genes-path', type=str, default='./data/genes_list.tsv')

    parser.add_argument('--expression-path', type=str, default='exp_array.tsv')
    parser.add_argument('--mutation-path', type=str, default='simple_somatic_mutation.open.tsv')
    parser.add_argument('--output-path', type=str, default='./output')

    run(parser.parse_args())
