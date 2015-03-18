# spectra-cluster

# Introduction
spectra-cluster is an open-source (Apache 2 licensed) library for clustering MS/MS spectra produced by mass spectrometers. It is designed with the goal of making
 spectra clustering easy to use and experiment for the community. It offers the following features out-of-box:

  - A collection of both classic and new algorithms for measuring [spectra similarities](https://github.com/spectra-cluster/spectra-cluster/tree/master/src/main/java/uk/ac/ebi/pride/spectracluster/similarity).
  - A set of [engines](https://github.com/spectra-cluster/spectra-cluster/tree/master/src/main/java/uk/ac/ebi/pride/spectracluster/engine) for clustering spectra together.
  - A set of [normalizers](https://github.com/spectra-cluster/spectra-cluster/tree/master/src/main/java/uk/ac/ebi/pride/spectracluster/normalizer) for normalising spectral peaks.
  - A set of [filters](https://github.com/spectra-cluster/spectra-cluster/tree/master/src/main/java/uk/ac/ebi/pride/spectracluster/util/predicate) and [functions](https://github.com/spectra-cluster/spectra-cluster/tree/master/src/main/java/uk/ac/ebi/pride/spectracluster/util/function) for pre-processing spectra, such as removing noisy peaks.
  - A set of cleanly defined data models and interfaces that represents spectra, peptide spectrum matches, and clusters.
  - Read in spectra and write out clustering results

# Getting started

### Installtion
You will need to have [Maven](http://maven.apache.org/) installed in order to build and use the spectra-cluster library.

Add the following snippets in your Maven pom file:

```maven
<!-- spectra-cluster dependency -->
<dependency>
    <groupId>uk.ac.ebi.pride.spectracluster</groupId>
    <artifactId>spectra-cluster</artifactId>
    <version>${current.version}</version>
</dependency>
```

```maven
 <!-- EBI repo -->
 <repository>
     <id>nexus-ebi-repo</id>
     <url>http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo</url>
 </repository>

 <!-- EBI SNAPSHOT repo -->
 <snapshotRepository>
    <id>nexus-ebi-repo-snapshots</id>
    <url>http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo-snapshots</url>
 </snapshotRepository>
```

### Running the library
TBD

# Getting help
If you have questions or need additional help, please contact the PRIDE Helpdesk at the EBI.

email: pride-support@ebi.ac.uk

# Feedback
Please give us your feedback, including error reports, suggestions on improvements, new feature requests. You can do so by opening a new issue at our [issues section](https://github.com/spectra-cluster/spectra-cluster/issues) 

# How to cite
Please cite this library using one of the following publications:
- Griss J, Foster JM, Hermjakob H, Vizcaíno JA. PRIDE Cluster: building the consensus of proteomics data. Nature methods. 2013;10(2):95-96. doi:10.1038/nmeth.2343. [PDF](http://www.nature.com/nmeth/journal/v10/n2/pdf/nmeth.2343.pdf),  [HTML](http://www.nature.com/nmeth/journal/v10/n2/full/nmeth.2343.html),  [PubMed](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3667236/)

# Contribute
TBD
# License
This project is available under the [Apache 2](http://www.apache.org/licenses/LICENSE-2.0.html) open source software (OSS) license.
