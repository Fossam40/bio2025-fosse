from collections import Counter
import matplotlib.pyplot as plt
import random

import mmh3

class BloomFilter:
    def __init__(self, size, nbhash):
        self.table = [0] * size
        self.nbhash = nbhash
        self.size = size
    
    def add(self, word):
        positions = [(mmh3.hash(word, i, signed = False) % self.size) for i in range(self.nbhash)]
        for i in positions:
            self.table[i] +=1
    
    def query(self, word):
        positions = [(mmh3.hash(word, i, signed = False) % self.size) for i in range(self.nbhash)]
        minseen = self.table[positions[0]]
        for i in positions:
            minseen = min(i, minseen)
        return minseen
    

def calculate_kmer_frequency(inputs, k, size_bloom, nbhash):
    counts = []
    bloom = BloomFilter(size_bloom, nbhash)
    for i in range(len(inputs)-k+1):
        kmer = inputs[i:i+k]
        nb_occ = bloom.query(kmer)
        bloom.add(kmer)
        if(len(counts)<=nb_occ):
            counts.append(1)
        else:
            counts[nb_occ] += 1
        
    return counts

def parse_fasta(filename):
    sequences = []
    with open(filename, "r") as file:
        for line in file:
            if(line[0] != ">"):
                sequences.append(line[:-1])
    return ''.join(sequences)

def main():
    sequence = parse_fasta("../github/Homework1/genome_hw1.fa")
    y = calculate_kmer_frequency(sequence, 10, 1000, 3)
    x = [i for i in range(len(y))]
    plt.figure(figsize=(8, 6))
    plt.plot(x, y)
    plt.title('Number of kmers agains number of appearances')
    plt.xlabel('Cumulative number of appearances')
    plt.ylabel('Number of kmers')
    #plt.show()
    plt.savefig("word_complexity.png")

if __name__ == "__main__":
    main()
