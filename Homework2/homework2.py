#! /usr/bin/env python3
from collections import Counter
import matplotlib.pyplot as plt


'''
For some reason, BLAST won't let me have results for sequences 1 and 4 even when I tweak the expected value.

Sequences are, in this order: , Sars-Cov-2, HIV1, , Sars-Cov-2 , Streptococcus (agalactiae, equi or suis)
'''

def parse_fasta(filename):
    sequences = []
    with open(filename, "r") as file:
        for line in file:
            if(line[0] != ">"):
                sequences.append(line[:-1])
    return ''.join(sequences)

letters_to_num = {'A': 1, 'C': 2, 'T': 4, 'G': 3}
letters_to_revnum = {'A': 4, 'C': 3, 'T': 1, 'G': 2}

def calculate_20mer_set(inputs):
    counts = set()
    current_read = 0
    current_rev_read = 0
    mod = (4**30)+1
    constant = 4**29
    for i in range(30):
        current_read = 4*current_read + letters_to_num[inputs[i]]
        current_rev_read = 4*current_rev_read + letters_to_revnum[inputs[29-i]]
    counts.add(min(current_read, current_rev_read))

    for i in range(31, len(inputs)):
        if(inputs[i]!="N" and inputs[i-20]!="N"):
            current_read = 4*(current_read % mod) + letters_to_num[inputs[i]]
            current_rev_read = current_rev_read - letters_to_revnum[inputs[i-20]] + letters_to_revnum[inputs[i]]*constant
            counts.add(min(current_read, current_rev_read))
    return counts

def jaccard_index(set1,  set2):
    return len(set1 & set2) / len(set1 | set2)

def kmer_counting(inputs):
    counts = [Counter() for i in range(29)]
    current_read = [0] * 29
    for i in range(len(inputs)):
        for j in range(min(30, i+1),2, -1):
            current_read[j-2] = current_read[j-3]*4 + letters_to_num[inputs[i]]
            counts[j-2][current_read[j-2]]+=1
        if(i>0):
            current_read[0] = 4*letters_to_num[inputs[i-1]] + letters_to_num[inputs[i]]
            counts[0][current_read[0]]+=1
    return [len(a) for a in counts]

def main():
    print("Computing Jaccard indexes")
    seqs = [""]*6

    for i in range(1,7):
        seqs[i-1] = parse_fasta("file" + str(i) + ".fa")

    kmerseqs = [calculate_20mer_set(a) for a in seqs]
    
    for i in range(1,7):
        for j in range(i+1, 7):
            print("     -> Index between sequences " + str(i) + " and " + str(j) + " : " + str(jaccard_index(kmerseqs[i-1], kmerseqs[j-1])))
    #All indexes are 0 so very similar kmer representation
    x = [i for i in range(2,31)]
    y = kmer_counting(seqs[2])
    plt.plot(x, y, marker='o', linestyle='-', label = "Real Genome")
    plt.title('K-mer count against K ')
    plt.xlabel('Value of K')
    plt.ylabel('Count of k-mer')
    plt.savefig("kmercount.png")

if __name__ == "__main__":
    main()
    
    


