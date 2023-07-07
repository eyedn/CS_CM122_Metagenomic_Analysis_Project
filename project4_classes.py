###############################################################################
#   Aydin Karatas
#   CS CM122 Project 4
#   project4_classes.py
###############################################################################
import project4_functions as f
from typing import Dict, List

# create a read object that knows it label and sequence
class Read:
    def __init__(self, label: str, sequence: str) -> None:
        self.lable = label
        self.sequence = sequence
        self.genome: Dict[Genome, int] = {}
        self.origin_genome: Genome = None

    def map_to_index(self, index: "Minimizer_Index", k: int, l: int) -> None:
        possible_genomes = {}
        for i in range(len(self.sequence) - k + 1):
            minimizer = f.minimize(self.sequence[i:i+k], l)
            # find which genomes are associated with this minimizer
            try:
                genome_list = index.index[minimizer]
            except:
                continue
            # record genomes this minimizer is found in
            for genome in genome_list:
                try:
                    possible_genomes[genome] += 1
                except:
                    possible_genomes[genome] = 1
        # assigne this read's genome as the one with the highest minimizer occ.
        self.genome = possible_genomes

    # given a list of which genomes should exists, assign a singe genome
    def identify_origin_genome(self, genomes: List["Genome"]) -> None:
        genome_in_list = False
        max_key = None
        max_value = float('-inf')
        for key, value in self.genome.items():
            if key in genomes:
                genome_in_list = True
                if value > max_value:
                    max_key = key
                    max_value = value
            elif not genome_in_list:
                if value > max_value:
                    max_key = key
                    max_value = value
        self.origin_genome = max_key

    def identify_genome_given_concensus_index(self) -> None:
        max_key = None
        max_value = float('-inf')
        for key, value in self.genome.items():
            if value > max_value:
                    max_key = key
                    max_value = value
        self.origin_genome = max_key



# create a genome object that knows it's label and sequence
class Genome:
    def __init__(self, label: str, sequence: str) -> None:
        self.label = label
        self.sequence = sequence

# create a minimzer index
class Minimizer_Index:
    def __init__(self) -> None:
        self.index: Dict[str, set(Genome)] = {}

    # add a genome to the index
    def add_genome(self, genome: "Genome", k: int, l: int) -> None:
        # iterate through kmers
        for i in range(len(genome.sequence) - k + 1):
            # find minimizer for each kmer
            kmer = genome.sequence[i:i+k]
            minimizer = f.minimize(kmer, l)
            # try to add the minizer to the index
            if minimizer in self.index.keys():
                if genome not in self.index[minimizer]:
                    self.index[minimizer].append(genome)
            # if failed, initalize indes with a list with this genome label 
            else: 
                self.index[minimizer] = [genome]