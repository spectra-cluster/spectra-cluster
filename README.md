# spectra-cluster
The core library for clustering MS spectra, including a number of similarity algorithms. 

# Main features
The main aim of the spectra-cluster library is to provide an open-source library for clustering 

# Similarity algorithms
Currently, we provide four algorithms for measuring the similarity between two MS2 spectra.
- [Dot product](https://github.com/spectra-cluster/spectra-cluster/blob/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity/FrankEtAlDotProduct.java) 
- [Hypergeometrics Test](https://github.com/spectra-cluster/spectra-cluster/blob/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity/HypergeometricScore.java)
- [Fisher Exact Test](https://github.com/spectra-cluster/spectra-cluster/blob/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity/FisherExactTest.java)
- [Intensity Ranking Correlation](https://github.com/spectra-cluster/spectra-cluster/blob/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity/IntensityRankCorrelation.java)

# Contribute
