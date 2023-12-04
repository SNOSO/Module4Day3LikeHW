import pandas as pd
import numpy as np
import random
from scipy.stats import percentileofscore

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

reverse_dict = {gene: key for key, value in FA_genes.items() for gene in value}

# random FA subnetwork with 1 gene per locus
def random_subnetworks(FA_genes):
    random_nets = {locus: random.choice(FA_genes[locus]) for locus in FA_genes.keys()}
    return [[random.choice(FA_genes[locus]) for locus in FA_genes.keys()] for _ in range(5000)]

prix_fixe = random_subnetworks(FA_genes)

def gene_connections(string, Gene1, Gene2):
    gene_connections = {}
    for gene1, gene2 in string[[Gene1, Gene2]].values.tolist():
        for gene in [gene1, gene2]:
            if gene1 in gene_connections: 
                if gene2 not in gene_connections[gene1]:
                    gene_connections[gene1].append(gene2) 
            else:
                gene_connections[gene1] = [gene2] 

            if gene2 in gene_connections: 
                if gene1 not in gene_connections[gene2]:
                    gene_connections[gene2].append(gene1) 
            else:
                gene_connections[gene2] = [gene1] 
    return gene_connections

gene_connections = gene_connections(string, 'Gene1', 'Gene2')

# Initialize a dictionary to store the counts
gene_counts = {}
for gene, connections in gene_connections.items():
    count = len(connections)
    gene_counts[gene] = count

# GENETIC ALGORITHM

# MUTATION STEP
def mutation(prix_fixe, FA_genes):
    mutated_subnetwork = []
    for subnetwork in prix_fixe:
        newsubnet = []
        for gene in subnetwork:
            if random.randint(1, 100) <= 5:
                locus = reverse_dict[gene]
                locus_size = len(FA_genes[locus])
                new_gene = FA_genes[locus][random.randint(0, locus_size - 1)]
                newsubnet.append(new_gene)
            else:
                newsubnet.append(gene)
        mutated_subnetwork.append(newsubnet)
    return mutated_subnetwork

mutated_subnetwork = mutation(prix_fixe, FA_genes)


# MATING STEP

class Subnetwork:
    def __init__(self, genes):
        self.genes = genes 
        self.edge_count = self.calculate_edge_count()
        self.selection_score = self.edge_count ** 3

    def calculate_edge_count(self):
        return sum(gene_counts.get(gene, 0) for gene in self.genes)

def normalize_scores(subnetworks):
    total_score = sum(sn.selection_score for sn in subnetworks)
    for sn in subnetworks:
        sn.selection_score /= total_score

def mate(subnetwork1, subnetwork2):
    genes = [random.choice((gene1, gene2)) for gene1, gene2 in zip(subnetwork1.genes, subnetwork2.genes)]
    return Subnetwork(genes)

def crossing(prix_fixe, FA_genes):
    subnetworks = [Subnetwork(subnetwork) for subnetwork in prix_fixe]
    normalize_scores(subnetworks)
    new_population = []
    for _ in range(5000):
        parents = random.choices(subnetworks, weights=[sn.selection_score for sn in subnetworks], k=2)
        new_population.append(mate(parents[0], parents[1]).genes)
    return new_population

child_subnetworks = crossing(prix_fixe, FA_genes)

def optimize(prix_fixe, child_subnetworks, FA_genes):
    prev_avg_density = sum(Subnetwork(subnetwork).edge_count for subnetwork in prix_fixe) / len(prix_fixe)
    while True:
        prix_fixe = mutation(prix_fixe, FA_genes)
        avg_density = sum(Subnetwork(subnetwork).edge_count for subnetwork in child_subnetworks) / len(child_subnetworks)
        if abs(avg_density - prev_avg_density) <= 0.5:
            break
        prev_avg_density = avg_density
    return prix_fixe

optimized_subnetworks = optimize(prix_fixe, child_subnetworks, FA_genes)

# Define the number of loci (L) using FA_genes keys and the number of genes in each locus (Gi) are the values in each key
L = len(FA_genes.keys())
Gi = [len(genes) for genes in FA_genes.values()]

# genes are the names of each value in each key
genes = [gene for genes in FA_genes.values() for gene in genes]

# Generate L random and disjoint sets of genes, with Gi genes per set i
gene_sets = []
for i in range(L):
    gene_sets.append(random.sample(genes, Gi[i]))

# use quantile based binning to bin genes into bins of equal size based on their counts
def bin_genes(gene_counts, bins):
    df = pd.DataFrame(list(gene_counts.items()), columns=['Gene', 'Connections'])
    bin_edges = range(0, df['Connections'].max() + bins, bins)
    labels = [f'{i}-{i+bins-1}' for i in bin_edges[:-1]]
    df['Bin'] = pd.cut(df['Connections'], bins=bin_edges, right=False, labels=labels)
    return df

bins = 10
bin_counts = bin_genes(gene_counts, bins)

# Replace each original candidate gene with a random gene selected from the same bin

# find FA genes in each bin
def find_genes(bin_counts, reverse_dict):
    mask = bin_counts['Gene'].isin(reverse_dict)
    fa_genes_in_bins = bin_counts.loc[mask, ['Gene', 'Bin']].set_index('Gene')['Bin'].to_dict()
    return fa_genes_in_bins

fa_genes_in_bins = find_genes(bin_counts, reverse_dict)

# find non-FA genes in each bin
def find_non_fa(bin_counts, reverse_dict):
    mask = bin_counts['Gene'].isin(reverse_dict)
    non_fa_genes_in_bins = bin_counts.loc[~mask, ['Gene', 'Bin']].set_index('Gene')['Bin'].to_dict()
    return non_fa_genes_in_bins

non_fa_genes_in_bins = find_non_fa(bin_counts, reverse_dict)

# define a function to replace elements in each prix_fixe subnetwork with genes in the same bins from non_fa_genes_in_bins
def replace_genes(prix_fixe, non_fa_genes_in_bins):
    return [set(random.choice(list(non_fa_genes_in_bins.keys())) for _ in range(len(prix_fixe[i]))) for i in range(len(prix_fixe))]

noninformative_network = replace_genes(prix_fixe, non_fa_genes_in_bins)

# Apply the genetic algorithm optimization method
noninformed_optimized_genes = optimize(noninformative_network, child_subnetworks, FA_genes)

# Calculate the average density of the noninformed_optimized_genes
def calculate_density(noninformed_optimized_genes):
    total_weight = 0
    max_possible_weight = 0
    for subnetwork in noninformed_optimized_genes:
        subnet_df = string.loc[(string['Gene1'].isin(subnetwork)) & (string['Gene2'].isin(subnetwork))]
        total_weight += subnet_df['Weight'].sum()  # Sum the weights of the edges in the subnetwork
        max_possible_weight += len(subnetwork) * (len(subnetwork) - 1) / 2
    return total_weight / max_possible_weight if max_possible_weight > 0 else 0

noninformed_average_density = calculate_density(noninformed_optimized_genes)
informed_average_density = calculate_density(optimized_subnetworks)

# compare test statistics for 1,000 random trials

def compare_statistics(net1, net2, trials):
    higher_count = 0
    for _ in range(trials):
        informed_density = calculate_density(net1)
        non_informed_density = calculate_density(net2)
        if non_informed_density > informed_density:
            higher_count += 1
    return higher_count / trials

p_value = compare_statistics(optimized_subnetworks, noninformed_optimized_genes, 1000)

# score genes in  based on their edge weights from FA_df
def gene_scoring(subnetwork, FA_genes, num_networks):
    def compute_density(subnetwork, full_network):
        subnet_df = full_network.loc[(full_network['Gene1'].isin(subnetwork)) & (full_network['Gene2'].isin(subnetwork))]
        return subnet_df['Weight'].sum()  # Sum the weights of the edges in the subnetwork

    gene_scores = {}
    loci = list(FA_genes.keys())
    full_densities = [compute_density(network, full_network) for network in subnetwork]
    prix_fixe = [set(network) for network in subnetwork]  # Convert to sets for faster membership tests

    for i in range(num_networks):
        if len(prix_fixe[i]) < len(loci):
            print(f"Network {i} has fewer genes than expected.")
            continue
        for j, locus in enumerate(loci):
            g_star = random.choice(list(subnetwork[i]))  # Select a random gene from the subnetwork
            new_network = subnetwork[i].copy() # make a new network with g_star removed
            new_network.remove(g_star) # remove g_star from the new network
            empty_density = compute_density(new_network, full_network) # calculate the density of the new network
            density_gi = full_densities[i] - empty_density  # calculate contribution of g* to the connectivity of the loci
            gene_scores[g_star] = density_gi / len(loci) # record contribution of g* to the connectivity of the loci
    return gene_scores

gene_scores = gene_scoring(subnetwork, full_network, 5000)

# convert gene_scores to a dataframe
gene_scores_df = pd.DataFrame.from_dict(gene_scores, orient='index', columns=['Score'])

# top 10 highest scored subnetworks 
def top_10_subnetworks(subnetwork, FA_df):
    def compute_density(subnetwork, full_network):
        subnet_df = full_network.loc[(full_network['Gene1'].isin(subnetwork)) & (full_network['Gene2'].isin(subnetwork))]
        return subnet_df['Weight'].sum()  # Sum the weights of the edges in the subnetwork

    full_densities = [compute_density(network, FA_df) for network in prix_fixe]
    prix_fixe = [set(network) for network in prix_fixe]  # Convert to sets for faster membership tests

    top_10 = []
    for i in range(10):
        top_density = max(full_densities)
        top_index = full_densities.index(top_density)
        top_network = prix_fixe[top_index]
        top_10.append((top_network, top_density))  # Append a tuple of the network and its score
        full_densities.remove(top_density)
        prix_fixe.remove(top_network)
    return top_10

top_10 = top_10_subnetworks(subnetwork, FA_df)

# convert top_10 to a dataframe
top_10_df = pd.DataFrame(top_10, columns=['Subnetwork', 'Score'])

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

# define a function to calculate p-values 
def calculate_p_values(gene_scores, num_permutations):
    p_values = {}
    all_scores = list(gene_scores.values())
    for gene, score in gene_scores.items():
        permuted_scores = [random.choice(all_scores) for _ in range(num_permutations)]
        p_value = percentileofscore(permuted_scores, score) / 100
        p_values[gene] = p_value
    return p_values

# calculate p-values
p_values = calculate_p_values(gene_scores, 1000)

# add p-values to the dataframe
final_ranking['P-value'] = pd.Series(p_values)

# calculate p-values for top_10 subnetworks
def calculate_p_values_top_10(top_10, num_permutations):
    p_values_top_10 = {}
    all_scores = [score for network, score in top_10]
    for network, score in top_10:
        permuted_scores = [random.choice(all_scores) for _ in range(num_permutations)]
        p_value = percentileofscore(permuted_scores, score) / 100
        p_values_top_10[tuple(network)] = p_value
    return p_values_top_10

pval_top_10 = calculate_p_values_top_10(top_10, 1000)

# create a file for each subnetwork in pval_top_10
def pval_output(pval_top_10):
    for i, (network, pval) in enumerate(pval_top_10.items()):
        with open(f'/Users/milligan-mcclellanlab/Desktop/Day3_Output_{i+1}_pval{pval}.txt', 'w') as f:
            for gene in network:
                f.write(f'{gene}\n')
    return

pval_output(pval_top_10)



