#! /usr/bin/env python3
from xopen import xopen
from readfa import readfq
import matplotlib
# We need, for every node, a list of outgoing edges (represented as the next letter). 
# We can store this in Python by using a dictionary with nodes as key and node lists as values.

def canonical_kmer(kmer):
    reverse_dict = {'A':'T', 'T': 'A', 'C':'G', 'G':'C'}
    reverse_kmer = ''.join([reverse_dict[e] for e in kmer[::-1]])
    return min(reverse_kmer, kmer)


def create_dbg(file, k, t):
    kmer_counts = {}
    total_kmers = 0
    kmers = set()
    with xopen(file) as fasta:
        for _,seq,_ in readfq(fasta):
            for i in range(len(seq)-k+1):
                kmer = canonical_kmer(seq[i:(i+k)])
                prefix = kmer[:-1]
                total_kmers+=1

                if(kmer in kmer_counts): 
                    kmer_counts[kmer]+=1
                else:
                    kmer_counts[kmer] = 1
                if(kmer_counts[kmer] >= t):
                    kmers.update([kmer])
    print(f"{total_kmers} kmers have been processed. {len(kmers)} distinct kmers are over the threshold.")
    return kmers


def unitig_from(dbg_graph, kmer):
    if(kmer not in dbg_graph):
        return ""
    visited = set([""])
    sufixes = []
    prefixes = []
    current_kmer = kmer
    while(current_kmer not in visited):
        visited.update([canonical_kmer(current_kmer)])
        possible_cont = ""
        for nucl in ["A", "C", "T", "G"]:
            possibility = "".join([current_kmer, nucl])[1:]
            if(canonical_kmer(possibility) in dbg_graph):
                if(possible_cont != ""):
                    break
                else:
                    possible_cont = possibility
        current_kmer = possible_cont
        if(len(possible_cont)>0):
            sufixes.append(current_kmer[-1])    
     
    current_kmer = kmer
    visited.remove(canonical_kmer(current_kmer))
    while(current_kmer not in visited):
        visited.update([canonical_kmer(current_kmer)])
        possible_cont = ""
        for nucl in ["A", "C", "T", "G"]:
            possibility = "".join([nucl, current_kmer])[:-1]
            if(canonical_kmer(possibility) in dbg_graph):
                if(possible_cont != ""):
                    break
                else:
                    possible_cont = possibility
        current_kmer = possible_cont
        if(len(possible_cont)>0):
            prefixes.append(possible_cont[-1])

    prefixes.reverse()
    prefixes = "".join(prefixes)
    sufixes = "".join(sufixes)

    return "".join([prefixes, kmer, sufixes])

if __name__ == "__main__":
    '''
    Q3: 
    153563 kmers are over the threshold, of which 153469 are distinct.

    Q4:
    We take k=100 so as to have the cleanest dBG possible.
    - perfect reads, forward strand results in 30000 kmers are over the threshold, of which 27264 are distinct (t = 1).
    - perfect reads, forward and reverse strands results in 30000 kmers are over the threshold, of which 27179 are distinct (t = 1).

    Q7:
    The length of the unitigs first increase and then decrease. For example, here are the results from reads/ecoli_sample_reads.fasta.gz

    Unitig len for t=1 : 82
    Unitig len for t=2 : 4202
    Unitig len for t=3 : 6292
    Unitig len for t=4 : 2243

    Because of read errors, at the beginning, the graph is very dense, so non-branching paths are short. By filtering,
    we reduce error-induced edges, therefore increasing the length of the unitig. However, by increasing too much the threshold, we end
    up removing too many edges, so the unitig length goes down.

    This seems to be supported by the fact that this phenomena is less present on perfect read files.
    '''

    #myg = create_dbg("ecoli_genome_150k.fa", 31, 1)
    filename = "reads/ecoli_sample_reads.fasta.gz"
    dbg_graph_t1 = create_dbg(filename,31,1)
    dbg_graph_t2 = create_dbg(filename,31,2)
    dbg_graph_t3 = create_dbg(filename,31,3)
    dbg_graph_t4 = create_dbg(filename,31,4)
    print("Results from " + filename)
    print("Unitig len for t=1 : " + str(len(unitig_from(dbg_graph_t1, "CGCTCTGTGTGACAAGCCGGAAACCGCCCAG"))))
    print("Unitig len for t=2 : " + str(len(unitig_from(dbg_graph_t2, "CGCTCTGTGTGACAAGCCGGAAACCGCCCAG"))))
    print("Unitig len for t=3 : " + str(len(unitig_from(dbg_graph_t3, "CGCTCTGTGTGACAAGCCGGAAACCGCCCAG"))))
    print("Unitig len for t=4 : " + str(len(unitig_from(dbg_graph_t4, "CGCTCTGTGTGACAAGCCGGAAACCGCCCAG"))))
            

    