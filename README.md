<div align="center">
  <img src="header.png"><br><br>
</div>
# Identification-of-brain-based-disorders

A MATLAB package for identifying neuropsychiatric disorders using unsupervised clustering on resting-state functional magnetic resonance imaging (rs-fMRI) connectivity features. Details are described in the paper "Investigating the correspondence of clinical diagnostic grouping with underlying neurobiological and phenotypic clusters using unsupervised learning". 

## Input data format

* clinical_diagnostic_labels: csv
* phenotypic_variables: csv
* connectivity_features: MAT-file with connectivity feature values `connectivity` and feature names `coordinates`.
```
# SFC: static functional connectivity
# DFC: variance of dynamic functional connectivity
# SEC: static effective connectivity
# DEC: variance of dynamic effective connectivity
e.g., 'DFC_19_64': DFC computed between 19 and 64 ROIs
```

## Clustering and feature selection

The main function for this package is defined in `runPipeline.m`, which takes a method name index, an input matrix, and an optional reference label.

### Input parameters

#### `method_name`: int

Indicate the selected clustering method. 0: hierarchical clustering, 1: DPC, and 2: OPTICS

#### `input_data`: matrix, shape = (n_samples, n_features)

#### `reference_label`: array, shape = (n_samples, 1), optional

If `reference_label` is `None`, use Calinski-Harabasz index as the criterion for 1) determination of optimal value of each input paramter of selected clustering method and 2) evaluation in feature selection. Otherwise, use cluster similarity between obtained label and reference label as the criterion for evaluation in feature selection.

### Returns

#### `bestScore`: float

#### `bestLabel`: array, shape = (n_samples, 1)

#### `bestFeature`: array, shape = (1, m), 0 < m < n_features

## Example

```
>> input_file = load('../data/connectivity_features/PTSD_connectivity.mat');
>> input_data = input_file.connectivities;

>> reference_label = load('sample_label.csv');

>> [bsetScore, bestLabel, bestFeature] = runPipeline(1, input_data, reference_label); % select DPC as clustering method
```

## References

[1] Rodriguez, A., Laio, A., 2014. Clustering by fast search and find of density peaks. Science 344, 1492-6. [Link](http://science.sciencemag.org/content/344/6191/1492.full)

[2] Liao, W., Chen, H., Yang, Q., Lei, X., 2008. Analysis of fMRI data using improved self-organizing mapping and spatio-temporal metric hierarchical clustering. IEEE Trans. Med. Imaging 27, 1472-1483. [Link](http://ieeexplore.ieee.org/abstract/document/4494444/)

[3] M. Ankerst, M. Breunig, H. Kriegel, J. Sander, 1999. OPTICS: Ordering Points To Identify the Clustering Structure. ACM SIGMOD Record, 28, 49-60. [Link] (http://dl.acm.org/citation.cfm?id=304187)

[4] Dy, J.G., Brodley, C.E., 2004. Feature Selection for Unsupervised Learning. J. Mach. Learn. Res. 5., 845-889. [Link] (http://www.ece.neu.edu/fac-ece/jdy/papers/dy04a.pdf)
