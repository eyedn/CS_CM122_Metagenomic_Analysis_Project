###############################################################################
#   Aydin Karatas
#   CS CM122 Project 4
#   project4_main.py
###############################################################################
import project4_functions as f
from sys import argv

def main():
    # my_dir = "sample"
    # output = "predictions_sample.csv"
    my_dir = argv[1]
    output = argv[2]

    # take in value for k
    input("Prove a k value. Press Enter to continue...")
    k = int(input("Enter the value for k: "))
    # take in value of l (minimizer size)
    input("Prove an l value for minimizer size. Press Enter to continue...")
    l = int(input("Enter the value for l: "))

    # read in genomes given the directory they exist in
    my_genomes, my_reads = f.read_genomes_and_reads(my_dir)
    
    # create minimizer index from the genomes
    minimizer_index = f.create_index(my_genomes, k, l)
    f.map_reads_to_index(my_reads, minimizer_index, k, l)
    f.find_concensus_genomes(my_reads)
    f.print_results(my_reads, output)
    
if __name__ == "__main__":
    main()