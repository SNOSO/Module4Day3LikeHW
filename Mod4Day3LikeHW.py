import pandas as pd
import numpy as np
import random

string = pd.read_csv('/Users/milligan-mcclellanlab/Desktop/STRING1.txt', delimiter='\t', names=['Gene1', 'Gene2', 'Weight'])

# define a function to create dictionary with locus as keys and genes as values
def fa_genes(file_name):
    FA_genes = {}
    for line in open('/Users/milligan-mcclellanlab/Desktop/Input.gmt.txt', 'r'):
        locus = line.split("\t")[0] 
        genes = line.split("\t")[2:]
        FA_genes[locus] = genes 
    return FA_genes

FA_genes = fa_genes('/Users/milligan-mcclellanlab/Desktop/Input.gmt.txt')

# define a function to create a network of FA genes
def FA_network(FA_genes, string, Gene1, Gene2, Weight):
    FA_genes = set(gene for genes in FA_genes.values() for gene in genes)
    mask = string[Gene1].isin(FA_genes) & string[Gene2].isin(FA_genes)
    network = string.loc[mask, [Gene1, Gene2, Weight]].values.tolist()
    return network

FA_network = FA_network(FA_genes, string, 'Gene1', 'Gene2', 'Weight')

# convert FA_network to a dataframe
FA_df = pd.DataFrame(FA_network, columns=['Gene1', 'Gene2', 'Weight'])

# random FA subnetwork with 1 gene per locus
def random_subnetworks(FA_genes):
    random_nets = {locus: random.choice(FA_genes[locus]) for locus in FA_genes.keys()}
    return [[random.choice(FA_genes[locus]) for locus in FA_genes.keys()] for _ in range(5000)]

prix_fixe = random_subnetworks(FA_genes)

# score genes in prix_fixe based on their edge weights from FA_df
def gene_scoring(prix_fixe, FA_genes, num_networks):
    def compute_density(subnetwork, full_network):
        subnet_df = full_network.loc[(full_network['Gene1'].isin(subnetwork)) & (full_network['Gene2'].isin(subnetwork))]
        return subnet_df['Weight'].sum()  # Sum the weights of the edges in the subnetwork

    gene_scores = {}
    loci = list(FA_genes.keys())
    full_densities = [compute_density(network, FA_df) for network in prix_fixe]
    prix_fixe = [set(network) for network in prix_fixe]  # Convert to sets for faster membership tests

    for i in range(num_networks):
        if len(prix_fixe[i]) < len(loci):
            print(f"Network {i} has fewer genes than expected.")
            continue
        for j, locus in enumerate(loci):
            g_star = random.choice(list(prix_fixe[i]))  # Select a random gene from the subnetwork
            new_network = prix_fixe[i].copy() # make a new network with g_star removed
            new_network.remove(g_star) # remove g_star from the new network
            empty_density = compute_density(new_network, FA_df) # calculate the density of the new network
            density_gi = full_densities[i] - empty_density  # calculate contribution of g* to the connectivity of the loci
            gene_scores[g_star] = density_gi / len(loci) # record contribution of g* to the connectivity of the loci
    return gene_scores

gene_scores = gene_scoring(prix_fixe, FA_genes, 5000)

# convert gene_scores to a dataframe
gene_scores_df = pd.DataFrame.from_dict(gene_scores, orient='index', columns=['Score'])

reverse_dict = {gene: key for key, value in FA_genes.items() for gene in value}

# final relative ranking of each gene in each locus
def final_ranking(gene_scores_df):
    gene_scores_df['Rank'] = gene_scores_df['Score'].rank(ascending=False)
    # add locus column
    gene_scores_df['Locus'] = gene_scores_df.index.map(reverse_dict)
    return gene_scores_df

final_ranking = final_ranking(gene_scores_df)

# use NA in Rank column if Score is 0
final_ranking['Rank'] = np.where(final_ranking['Score'] == 0, np.nan, final_ranking['Rank'])

# export final_ranking to a txt file
final_ranking.to_csv('/Users/milligan-mcclellanlab/Desktop/Day3_Output.gmt', sep='\t', index=False)

 