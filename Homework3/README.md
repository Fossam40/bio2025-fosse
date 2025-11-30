# Read assembly
This repository has functions to create a de Bruijn graph, and the unitigs associated to it, from a set of reads. Here are the answers to this problem set's questions:
    
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
    
