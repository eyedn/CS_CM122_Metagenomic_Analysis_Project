###############################################################################
#   Aydin Karatas
#   CS CM122 Project 4
#   project4_functions.py
###############################################################################
import project4_classes as c
import datetime
from typing import List, Dict
import os
import math

# update variables for printing to terminal
update_file = 100
update_genome = 100
update_reads = 10000

# given a kmer, find its minimizer:
def minimize(kmer: str, l: int) -> None:
    minimizer = None
    for j in range(len(kmer) - l + 1): 
        if minimizer == None or minimizer > kmer[j:j+l]:
            minimizer = kmer[j:j+l]
    return minimizer

# given a dictionary [str : int], return the str with the high int
def get_key_with_max_value(dictionary: Dict[c.Genome, int], genomes: List[c.Genome]) -> c.Genome:
    genome_in_list = False
    max_key = None
    max_value = float('-inf')
    for key, value in dictionary.items():
        if key in genomes:
            genome_in_list = True
            if value > max_value:
                max_key = key
                max_value = value
        elif not genome_in_list:
            max_key = key
            max_value = value
    return max_key

# read in all genomes and create minimizer index:
def read_genomes_and_reads(directory: str) -> None:
    print(f"{datetime.datetime.now()}: reading in data")
    genomes = []
    reads = []
    update = 0
    # Iterate over each file in the directory
    for filename in os.listdir(directory):
        if update != 0 and update % update_file == 0:
            print(f"{datetime.datetime.now()}: {update} files read in")
        # Get the absolute file path
        file_path = os.path.join(directory, filename)
        # Check if the path corresponds to a fasta file
        if os.path.isfile(file_path) and filename.endswith('.fasta'):
            # if we encounter the reads, file, we will read in the reads
            if filename.endswith('reads.fasta'):
                with open(file_path, 'r') as r:
                    for line in r.readlines():
                        if line.startswith(">"):
                            reads.append(c.Read(line.strip(), ""))
                            continue
                        reads[-1].sequence += line.strip()
            # if we encounter a genome file, we will read in the gneome
            else:
                with open(file_path, 'r') as g:
                    for line in g.readlines():
                        if line.startswith(">"):
                            genomes.append(c.Genome(line.strip(), ""))
                            continue
                        genomes[-1].sequence += line.strip()
            update += 1
    print(f"{datetime.datetime.now()}: genomes and reads have been read in")
    return genomes, reads

# create minimizer index given all the genomes:
def create_index(genomes: List[c.Genome], k: int, l: int) -> c.Minimizer_Index:
    print(f"{datetime.datetime.now()}: creating index")
    index = c.Minimizer_Index()
    for i, genome in enumerate(genomes):
        if i != 0 and i % update_genome == 0:
            print(f"{datetime.datetime.now()}: {i} of {len(genomes)} indexed")
        index.add_genome(genome, k, l)
    print(f"{datetime.datetime.now()}: indexing complete")
    return index

# sample 20% of reads to find concensus genomes and adjust index
def sample_reads(reads: List[c.Read], index: c.Minimizer_Index, k: int, l: int, 
                 sample_size: float) -> c.Minimizer_Index:
    print(f"{datetime.datetime.now()}: sampling first {sample_size} of reads")
    last_index = int(math.floor(len(reads)*sample_size))
    for i, read in enumerate(reads):
        if i == last_index:
            break
        if i != 0 and i % update_reads == 0:
            print(f"{datetime.datetime.now()}: {i} of {len(reads)} reads sampled")
        read.map_to_index(index, k, l)
    print(f"{datetime.datetime.now()}: finding concensus genomes")
    # count genome occurances
    genome_counts = {}
    for i, read in enumerate(reads):
        if i == last_index:
            break
        for genome, count in read.genome.items():
            if genome in genome_counts:
                genome_counts[genome] += count
            else:
                genome_counts[genome] = count
    sorted_genomes = sorted(genome_counts.items(), key=lambda x: x[1], reverse=True)
    # print distribution of genomes to the user
    with open("distribution.txt", 'w') as file:
        for genome, count in sorted_genomes:
            file.write(f"Genome: {genome.label}\tCount: {count}\n")
    # the user will pick how many genomes to consider
    input("Prove how mant genomes will be considered. Press Enter to continue...")
    k = int(input("Enter the value the number of genomes to be considered: "))
    top_genomes = [genome for genome, _ in sorted_genomes[:k]]
    print(f"{datetime.datetime.now()}: adjusting index to only contain top genomes")
    # adjust index to only include top genomes
    keys_to_delete = []
    for minimizer in list(index.index):
        intersect = set(top_genomes) & set(index.index[minimizer])
        if len(intersect) == 0:
            keys_to_delete.append(minimizer)
    for key in keys_to_delete:
        del index.index[key]
    print(f"{datetime.datetime.now()}: mapping complete")
    return index

# map all reads in a minimzer index:
def map_reads_to_index(reads: List[c.Read], index: c.Minimizer_Index, k: int, 
                       l: int) -> None:
    print(f"{datetime.datetime.now()}: mapping reads to index")
    for i, read in enumerate(reads):
        if i != 0 and i % update_reads == 0:
            print(f"{datetime.datetime.now()}: {i} of {len(reads)} reads mapped")
        read.map_to_index(index, k, l)
    print(f"{datetime.datetime.now()}: mapping complete")

# after mapping reads, find which genomes are most commonly found
def find_concensus_genomes(reads: List[c.Read]) -> None:
    print(f"{datetime.datetime.now()}: finding concensus genomes")
    genome_counts = {}
    for read in reads:
        for genome, count in read.genome.items():
            if genome in genome_counts:
                genome_counts[genome] += count
            else:
                genome_counts[genome] = count
    sorted_genomes = sorted(genome_counts.items(), key=lambda x: x[1], reverse=True)
    # print distribution of genomes to the user
    with open("distribution.txt", 'w') as file:
        for genome, count in sorted_genomes:
            file.write(f"Genome: {genome.label}\tCount: {count}\n")
    # the user will pick how many genomes to consider
    input("Prove how mant genomes will be considered. Press Enter to continue...")
    k = int(input("Enter the value the number of genomes to be considered: "))
    top_genomes = [genome for genome, _ in sorted_genomes[:k]]
    print(f"{datetime.datetime.now()}: consensus genomes found")
    print(f"{datetime.datetime.now()}: identifying correct genomes for each read")
    # for all reads, idenitfy their tru genomes
    for i, read in enumerate(reads):
        if i != 0 and i % update_reads == 0:
            print(f"{datetime.datetime.now()}: {i} of {len(reads)} read mapped to one genome")
        read.identify_origin_genome(top_genomes)
    print(f"{datetime.datetime.now()}: correct genomes for each read")

# map all reads to in a minimzer index:
def assign_genomes(reads: List[c.Read], index: c.Minimizer_Index, k: int, 
                       l: int) -> None:
    print(f"{datetime.datetime.now()}: mapping reads to index")
    for i, read in enumerate(reads):
        if i != 0 and i % update_reads == 0:
            print(f"{datetime.datetime.now()}: {i} of {len(reads)} reads mapped")
        read.map_to_index(index, k, l)
        read.identify_genome_given_concensus_index()
    print(f"{datetime.datetime.now()}: mapping complete")

# print results of which genome a read came from
def print_results(reads: List[c.Read], output: str) -> None:
    print(f"{datetime.datetime.now()}: generating results")
    with open(output, "w") as f:
        for i, read in enumerate(reads):
            if i != 0 and i % update_reads == 0:
                print(f"{datetime.datetime.now()}: {i} of {len(reads)} results processed")
            if read.origin_genome == None:
                f.write(f"{read.lable}\tNONE\n")
            else:
                genome_label = read.origin_genome.label[1:].split(' ')[0]
                f.write(f"{read.lable}\t{genome_label}\n")
    print(f"{datetime.datetime.now()}: finished!")