from collections import Counter
import matplotlib.pyplot as plt
import random

letters_to_num = {'A': 1, 'C': 2, 'T': 3, 'G': 4}

# Takes a read as argument, returns a 29 long array with the number of k-mers, k<=30
def calculate_word_complexity(inputs):
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

# Very simple fasta parsing
def parse_fasta(filename):
    sequences = []
    with open(filename, "r") as file:
        for line in file:
            if(line[0] != ">"):
                sequences.append(line[:-1])
    return ''.join(sequences)

# Generates Fibonacci words
def generate_fibonacci(length):
    init = ['A']
    while(len(init)<length):
        new_init = []
        for i in init:
            new_init.append('A')
            if(i=='A'):
                new_init.append('T')
        init = new_init
    return ''.join(init)
    
#Generates random words
def generate_random(length):
    seq = [random.choice(['A','C','G','T']) for i in range(length)]
    return ''.join(seq)

# Example plots, comparing Fibonacci vs Random vs actual sequences
def main():
    sequence = parse_fasta("genome_hw1.fa")
    x = [i for i in range(2,31)]
    y0 = calculate_word_complexity(sequence)
    y1 = calculate_word_complexity(generate_fibonacci(len(sequence)))
    y2 = calculate_word_complexity(generate_random(len(sequence)))
    plt.figure(figsize=(8, 6))
    plt.plot(x, y0, marker='o', linestyle='-', label = "Real Genome")
    plt.plot(x, y1, marker='o', linestyle='-', label = "Fibonacci")
    plt.plot(x, y2, marker='o', linestyle='-', label = "Random")
    plt.title('Subword complexity against subword size')
    plt.xlabel('Subword length K')
    plt.ylabel('F(k)')
    plt.legend()
    plt.savefig("word_complexity.png")

if __name__ == "__main__":
    main()
