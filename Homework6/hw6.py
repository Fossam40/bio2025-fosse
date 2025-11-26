#! /usr/bin/env python3
from xopen import xopen
from readfa import readfq
import matplotlib
import numpy as np
import csv

def reverse_read(read):
    reverse_dict = {'A':'T', 'T': 'A', 'C':'G', 'G':'C'}
    reverse_read = ''.join([reverse_dict[e] for e in read[::-1]])
    return reverse_read

def extract_kmers(seq, length):
    kmer_to_pos = {}
    current_kmer = 0
    mult_factor = 4**(length - 1)
    letter_to_num = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    for i in range(length):
        current_kmer *= 4
        current_kmer += letter_to_num[seq[i]]
    kmer_to_pos[current_kmer] = {0}

    for i in range(length, len(seq)):
        current_kmer -= mult_factor*letter_to_num[seq[i - length]]
        current_kmer *= 4
        current_kmer += letter_to_num[seq[i]]
        if current_kmer in kmer_to_pos:
            kmer_to_pos[current_kmer].add(i)
        else:
            kmer_to_pos[current_kmer] = {i}
    return kmer_to_pos

def edDistDp(x, y):
    """ Calculate edit distance between sequences x and y using
    matrix dynamic programming. Return distance. """
    dp = np.zeros((len(x)+1, len(y)+1), dtype=int)
    dp[0, 1:] = range(1, len(y)+1)
    dp[1:, 0] = range(1, len(x)+1)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            dp[i, j] = min(dp[i-1, j-1]+delt, dp[i-1, j]+1, dp[i, j-1]+1)
    return dp[len(x), len(y)]

# Takes a read, the beginning of a position, a reference and calculates 
# the edit distance between a portion of the ref and the read
def positionalEditDist(read, pos_beg, ref):
    subseq = ref[pos_beg: min(len(ref)-1, pos_beg + len(read) + len(read)//2)]
    return edDistDp(read, subseq)


def query_read(read, ref_kmer, ref, length):
    current_kmer = 0
    mult_factor = 4**(length - 1)
    letter_to_num = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

    for i in range(length):
        current_kmer *= 4
        current_kmer += letter_to_num[read[i]]

    min_dist = float("inf")
    min_pos = -1
    if(current_kmer in ref_kmer):
        for match_pos in ref_kmer[current_kmer]:
            pos_beg = max(0, match_pos - i)
            edist = positionalEditDist(read, pos_beg, ref)
            if(edist < min_dist):
                min_dist = edist
                min_pos = pos_beg

    seen_positions = set()

    for i in range(length, len(read)):
        current_kmer -= mult_factor*letter_to_num[ref[i - length]]
        current_kmer *= 4
        current_kmer += letter_to_num[ref[i]]
        if(current_kmer in ref_kmer):
            for match_pos in ref_kmer[current_kmer]:
                pos_beg = max(0, match_pos - i - len(read)//2)
                if(pos_beg not in seen_positions):
                    seen_positions.add(pos_beg)
                    edist = positionalEditDist(read, pos_beg, ref)
                    if(edist < min_dist):
                        min_dist = edist
                        min_pos = pos_beg
    return min_pos, min_dist

reads = "reads.fq.gz"
refs = 'ref.fa.gz'
k = 9

with open('output.tsv', 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    writer.writerow(["Query", "RefName", "Position", "Strand", "Score"])
    with xopen(refs) as fastarefs:
        for refname,ref,_ in readfq(fastarefs):
            ref_kmer = extract_kmers(ref, k)
            with xopen(reads) as fastareads:
                for readname,read,_ in readfq(fastareads):
                    forward = query_read(read, ref_kmer, ref, k)
                    reverse = query_read(reverse_read(read), ref_kmer, ref, k)
                    if(forward[1] < reverse[1]):
                        writer.writerow([readname, refname, forward[0], "+", forward[1]])
                    elif(forward[1] > reverse[1]):
                        writer.writerow([readname, refname, reverse[0], "-", reverse[1]])





        
