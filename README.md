#### CS CM 122 Project 4
### Aydin Karatas
I created a metagenomic tool to identify what genomes are present in a sample from a list of possible genomes. There were 5000 total genome in the list.

#### Setup:
We have a directory contain genome files in the form <*genome*.fasta>. 
This directory also contains a reads file in the form <*reads.fasta>.

#### Command:
To find which genome each read comes from, run the command...
    python3 project4a_main.py <directory_to_genomes_and_reads> <predictions.csv>

#### Program Overview:
The program will first create an index of minimizers in the from of a hashtable. The keys are minimizers based on kmers and the values are which genomes this minimizer was found in. The program will first prompt the user to provide a <l> value for the size of the sliding window and an <l> value that will be the minimizer size value. For this project, I chose <k=20> and <l=10>.

The program will read in all files first. Then, it will generate the index, which is the longest step. Once the index is genreated, the reads will be split up into minimizers based on kmers. The program will produce a file called <distribution.txt> that provides how many times a given genome was found in any of the minimizers of any read. This will prompt the user to provide how many genomes to include in the final query. For this project, I chose to only consider the top 7 most expressed genomes.

#### Output:
The output file for this code was provided as the second argument to the inital command. In this case, it was <predictions.csv>.