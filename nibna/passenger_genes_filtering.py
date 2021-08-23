from collections import defaultdict
import logging
from os.path import join
from typing import Dict, List, Set

# Libraries for matrix computations
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Libraries for passenger genes filtering
from pygam import LinearGAM, LogisticGAM, s, f, te
import random
import math

DATA_DIR = "C:\\Users\\Hoang Tung\\Documents\\GitHub\\NIBNAPlus\\Data\\"
OUT_DIR = "C:\\Users\\Hoang Tung\\Documents\\GitHub\\NIBNAPlus\\Data\\Output\\"

def CCDS_data_processing(id_file_name: str, length_file_name: str):
    # read in the CCDS current file
    CCDS_current = pd.read_csv(join(DATA_DIR, id_file_name), sep="\t")
    CCDS_current = CCDS_current[['gene','gene_id', 'ccds_id']]

    # process the CCDS protein data
    ccds_length = dict()
    with open(join(DATA_DIR,length_file_name), "r") as a_file:
        for line in a_file:
            if line[0] == '>':
                index = line.find('|')
                ccds_id = line[1:index]
                ccds_length[ccds_id] = 0
            else:
                aa_partial_length = len(line) - 1
                ccds_length[ccds_id] = ccds_length[ccds_id] + 3*aa_partial_length

    data_items = ccds_length.items()
    data_list = list(data_items)
    CCDS_cDNA_length = pd.DataFrame(data_list)

    # joint two table based in ccds_id
    CCDS_data = pd.merge(CCDS_current, CCDS_cDNA_length, how="left", left_on=['ccds_id'], right_on = [0]).drop(columns=[0]).rename(columns={1:"cDNA length"})
    return CCDS_data

def TCGA_data_processing(file_name: str, CCDS_data: pd.DataFrame):
    #read TCGA dataset
    TCGA = pd.read_csv(join(DATA_DIR,file_name), sep="\t")
    TCGA = pd.merge(TCGA, CCDS_data, how="left", left_on=['Entrez_Gene_Id'], right_on = ['gene_id'])
    return TCGA

def sample_filter(TCGA, CCDS_data, sample_id):
    # generate the mutation status for every CCDS gene for each sample
    mutation_data = TCGA.loc[TCGA['Tumor_Sample_Barcode'] == sample_id].reset_index()
    mutation_list = mutation_data['Entrez_Gene_Id'].tolist()
    all_ccds_gene = CCDS_data.copy()
    all_ccds_gene['mutation status'] = 0
    for i in range(len(all_ccds_gene)):
        gene_id = all_ccds_gene.loc[i]['gene_id']
        if gene_id in mutation_list:
            all_ccds_gene.at[i, 'mutation status'] = 1

    # prepare the dataset
    all_ccds_gene = all_ccds_gene.dropna() #drop genes without reported cDNA length
    all_ccds_gene = all_ccds_gene.sort_values(by='cDNA length', ascending=False)
    all_ccds_gene = all_ccds_gene.drop_duplicates(subset='gene', keep="first")

    gene_X = all_ccds_gene['cDNA length']
    gene_Y = all_ccds_gene['mutation status']

    # Fit using logistic GAM. F is monotonic increase cubic spline with 6 knots
    lgam = LogisticGAM(s(0, n_splines = 7, spline_order = 4, constraints = 'monotonic_inc')).fit(gene_X, gene_Y)

    # Find the PMV with reverse logit function
    fitted_value = lgam.partial_dependence(term=0, X=gene_X)
    pmv = [math.exp(x)/(1+math.exp(x)) for x in fitted_value]
    all_ccds_gene['PMV'] = pmv

    # resample
    all_ccds_gene['Number of sample contain'] = [0] *len(all_ccds_gene) #store the number of sample containing that gene as mutatated gene
    gene_list = all_ccds_gene['gene'].tolist()

    for i in range(1000):
        sample_gene = all_ccds_gene.sample(n=max([50, len(mutation_list)]), replace=False, weights=pmv)
        selected_genes = sample_gene['gene'].tolist()
        for gene in selected_genes:
            all_ccds_gene.loc[all_ccds_gene['gene'] == gene, 'Number of sample contain'] += 1

    return all_ccds_gene

def passenger_genes_filtering(TCGA: pd.DataFrame, CCDS_data: pd.DataFrame, CCG_genes:  pd.Series, refilter):
    """
    This function takes in the mutation status from TCGA breast cancer dataset and gene length from CCDS.
    It will perform the logistic general additive model to obtain mutation probability based on gene length
    for each genes. Then it will resample the dataset 1000 times to eliminate passenger mutation (incident > 50).
    Args:
        TCGA (pd.DataFrame): mutation status from TCGA breast cancer dataset
        CCDS_data (pd.DataFrame): gene length from CCDS. Currently, it is p38 version.
        CCG (pd.DataFrame): Cancer Census Gene list to cross validate
    """
    if refilter == True:
        sample_list = set(TCGA['Tumor_Sample_Barcode'].tolist())
        passenger_gene = list()
        for sample_id in sample_list:
            all_ccds_gene = sample_filter(TCGA, CCDS_data, sample_id)
            final_mut_gene = all_ccds_gene.loc[all_ccds_gene['mutation status'] == 1].loc[all_ccds_gene['Number of sample contain'] >= 50]
            filtered_gene = final_mut_gene['gene'].tolist()
            passenger_gene = passenger_gene + filtered_gene

        passenger_gene = list(set(passenger_gene))
        CCG_genes = CCG_genes.tolist()
        passenger_gene = np.setdiff1d(passenger_gene,CCG_genes)
    elif refilter == False:
        pgfile = pd.read_csv(join(DATA_DIR, 'passenger_gene.csv'))
        passenger_gene = pgfile['gene'].tolist()
        CCG_genes = CCG_genes.tolist()
        passenger_gene = np.setdiff1d(passenger_gene,CCG_genes)
    return passenger_gene

def critical_nodes_filter(critical_nodes, passenger_genes):
    original_length = len(critical_nodes)
    critical_nodes = critical_nodes[~critical_nodes['node'].isin(passenger_genes)].reset_index()
    diff = original_length - len(critical_nodes)
    print("Filtered " + str(diff) + " passenger genes from data")
    return critical_nodes
