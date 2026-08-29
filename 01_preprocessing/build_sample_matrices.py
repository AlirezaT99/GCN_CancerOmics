import os
import argparse

import numpy as np
import pandas as pd

from tqdm import tqdm

tqdm.pandas()

DATA_DIR = './data'
OUTPUT_DIR = './output'

cancers_with_ILMN = ['Pancreas']  # gene_id: ILMN_
cancers_with_NM = ['Nervous System']  # gene_id: NM_/NR/_
cancers_with_ENS_version = ['Blood']


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


def get_expression_data(cancer_type, expression_path):
    columns = ['icgc_donor_id', 'gene_id']
    data = read_chunk_by_chunk(cancer_type, expression_path, columns)
    return data.dropna()


def get_genes(genes_path, gene_class=None):
    genes = pd.read_csv(genes_path, sep='\t', header=0)
    genes.gene_symbol = list(map(lambda g: g[1:-1], genes.gene_symbol))
    if gene_class:
        genes = genes[genes['gene_class'] == gene_class]
    genes = genes[['gene_name', 'gene_symbol']] \
        .rename({'gene_name': 'gene_ensembl_id'}, axis=1)
    return genes


def store_summary(df, output_path, file_name):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    save_tsv(df, OUTPUT_DIR, f'result-{file_name}.tsv')
    print(f'>>> Summary for {file_name} (considering mutation and expression):')
    print('\tDonors in Common:', df.shape[0])
    print('\tGenes in Common:', df.shape[1])


def perform_analysis(args, cancer_type):
    genes = get_genes(args.genes_path)

    ### Mutation
    mut = get_mutation_data(cancer_type, args.mutation_path, mutation_type='single base substitution')
    mut_data = mut.rename({'gene_affected': 'gene_ensembl_id'}, axis=1)
    sign_mut_samples = pd.merge(genes, mut_data, how='left', on='gene_ensembl_id') \
        .drop(columns=['gene_ensembl_id']) \
        .drop_duplicates() \
        .dropna()

    ### Expression
    expr = get_expression_data(cancer_type, args.expression_path)

    #### Before this part the R script needs to have been run to convert Illumina probe to gene
    if cancer_type in cancers_with_ILMN:
        ILMN_genes = pd.read_csv(f'./{DATA_DIR}/{cancer_type}/converted_genes.csv')['Gene']
        expr['gene_symbol'] = ILMN_genes
        expr = expr.dropna()
    if cancer_type in cancers_with_ENS_version + cancers_with_NM:
        converted_genes = pd.read_csv(f'./{DATA_DIR}/{cancer_type}/converted_genes.csv')
        expr = pd.merge(expr, converted_genes, how='left', left_on="gene_id", right_on='initial_id')
        expr = expr.rename({'Gene': 'gene_symbol'}, axis=1)
        expr = expr.dropna()
    else:
        expr = expr.rename({'gene_id': 'gene_symbol'}, axis=1)

    ## Merge datasets
    ### Find intersection
    cols = ['icgc_donor_id', 'gene_symbol']
    mut_data = sign_mut_samples[cols]
    expr_data = expr[cols]

    common_donors = np.intersect1d(mut_data[['icgc_donor_id']], expr_data[['icgc_donor_id']])
    print('Initial common donors:', common_donors.shape)
    common_genes = np.intersect1d(mut_data[['gene_symbol']], expr_data[['gene_symbol']])
    print('Initial common genes:', common_genes.shape)

    # Narrow down both datasets
    final_mut = pd.merge(pd.Series(common_genes, name='gene_symbol'), mut_data, how='left', on='gene_symbol')
    final_mut = pd.merge(pd.Series(common_donors, name='icgc_donor_id'), final_mut, how='left', on='icgc_donor_id')
    final_expr = pd.merge(pd.Series(common_genes, name='gene_symbol'), expr_data, how='left', on='gene_symbol')
    final_expr = pd.merge(pd.Series(common_donors, name='icgc_donor_id'), final_expr, how='left', on='icgc_donor_id')

    updated_common_genes = np.intersect1d(final_expr.gene_symbol.unique(), final_mut.gene_symbol.unique())
    updated_common_donor = np.intersect1d(final_expr.icgc_donor_id.unique(), final_mut.icgc_donor_id.unique())
    final_mut = pd.merge(pd.Series(updated_common_genes, name='gene_symbol'), final_mut, how='left', on='gene_symbol')
    final_mut = pd.merge(pd.Series(updated_common_donor, name='icgc_donor_id'), final_mut, how='left',
                         on='icgc_donor_id')
    final_expr = pd.merge(pd.Series(updated_common_genes, name='gene_symbol'), final_expr, how='left', on='gene_symbol')
    final_expr = pd.merge(pd.Series(updated_common_donor, name='icgc_donor_id'), final_expr, how='left',
                          on='icgc_donor_id')

    final_mut.sort_values(by='gene_symbol', inplace=True)
    final_expr.sort_values(by='gene_symbol', inplace=True)
    sorted_common_genes = np.sort(updated_common_genes)

    ### Matrix generation
    result_mut = pd.DataFrame(index=updated_common_donor, columns=updated_common_genes).fillna(0)
    for idx, row in tqdm(final_mut.drop_duplicates().groupby('icgc_donor_id')['gene_symbol'].apply(list).iteritems()):
        result_mut.loc[idx, row] = 1

    result_expr = pd.DataFrame(index=updated_common_donor, columns=updated_common_genes).fillna(0)
    for idx, row in tqdm(final_expr.drop_duplicates().groupby('icgc_donor_id')['gene_symbol'].apply(list).iteritems()):
        result_expr.loc[idx, row] = 1

    result_values = result_expr.values * result_mut.values
    result = pd.DataFrame(data=result_values, index=updated_common_donor, columns=updated_common_genes)

    #### Store results
    store_summary(result, OUTPUT_DIR, f'result-{cancer_type}.csv')


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

    parser.add_argument('--cancer-type', type=str, default='Test')
    parser.add_argument('--run-all', type=bool, default=False)

    parser.add_argument('--data-path', type=str, default='./data')
    parser.add_argument('--genes-path', type=str, default='./data/genes_list.tsv')

    parser.add_argument('--expression-path', type=str, default='exp_array.tsv')
    parser.add_argument('--mutation-path', type=str, default='simple_somatic_mutation.open.tsv')
    parser.add_argument('--output-path', type=str, default='./output')

    run(parser.parse_args())
