# spectra-cluster

# Introduction
The core library for clustering MS spectra, including a number of similarity algorithms. 

# Features
The main aim of the spectra-cluster library is to provide an open-source library for clustering 

## Similarity algorithms

Currently, we provide four algorithms for measuring the similarity between two MS2 spectra.
- [Dot product](https://github.com/spectra-cluster/spectra-cluster/blob/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity/FrankEtAlDotProduct.java) 
- [Hypergeometrics Test](https://github.com/spectra-cluster/spectra-cluster/blob/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity/HypergeometricScore.java)
- [Fisher Exact Test](https://github.com/spectra-cluster/spectra-cluster/blob/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity/FisherExactTest.java)
- [Intensity Ranking Correlation](https://github.com/spectra-cluster/spectra-cluster/blob/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity/IntensityRankCorrelation.java)

# Getting started

### Getting spectra-cluster library
### How to use spectra-cluster library

# Getting help
If you have questions or need additional help, please contact the PRIDE Helpdesk at the EBI: pride-support at ebi.ac.uk (replace at with @).

# Giving your feedback
Please give us your feedback, including error reports, suggestions on improvements, new feature requests. You can do so by opening a new issue at our [issues section](https://github.com/spectra-cluster/spectra-cluster/issues) 

# License
This project is available under the [Apache 2](http://www.apache.org/licenses/LICENSE-2.0.html) open source software (OSS) license.

# How to cite
Please cite this library using one of the following publications:
- Griss J, Foster JM, Hermjakob H, Vizcaíno JA. PRIDE Cluster: building the consensus of proteomics data. Nature methods. 2013;10(2):95-96. doi:10.1038/nmeth.2343. [PDF](http://www.nature.com/nmeth/journal/v10/n2/pdf/nmeth.2343.pdf),  [HTML](http://www.nature.com/nmeth/journal/v10/n2/full/nmeth.2343.html),  [PubMed](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3667236/)

# Contribute
We welcome any contribution to the project
