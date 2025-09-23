# Exercice A
- For high quality reads (short and long), dBG are more suited, as the asymptotic complexity is better than OGs, and the lack of errors will not create many twigs and bubbles.

- For low quality reads (short), OGs are better as approximate overlaps can be used to partially compensate the errors. For long reads, dBG graphs with cleaning steps can be used to obtain somewhat correct results while avoiding a quadratic complexity.


# Exercice B

- We start with a multiset (a set where elements have a cardinality) of k-mer reads. We assume that the number of read errors is small.

- Build a dBG: For each of the kmers, compute a hash of the prefix and suffix (k-1)-mers. Then, for each suffix kmer, add an edge between the k-mers that have that kmer as suffix and the ones that have it as prefix. Complexity: O(k*(#reads)^2*log(#reads)) if we use a set for lookup.

- Burst dbG bubbles and twigs: Detect them by finding connex components. If a node has two or more outgoing edges and their connex components are different, there is a bubble or a twig. We keep the path whose k-mers are more present, ie the sum of the occurrences of its kmers is higher. We can precalculate this number by considering the DAG graph of the connex components and use dynamic programming over it, giving a total complexity of O(#(distinct reads))

- Find all unitigs in the resulting graph and output them: O(#(distinct reads))