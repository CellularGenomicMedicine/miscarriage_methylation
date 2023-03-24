# miscarriage_methylation
Methylation analysis of EM and CV samples from first trimester pregnancy loss to conduct cell composition deconvolution and to determine whether the sampled tissue types are distinct and whether they correspond to the expected cell sample types.

## Pre-processing & sample QC
1. [SWAN normalisation and sample / site QC](scripts/preprocessingSWAN.R)
    (the SNP heatmap can be extracted from the generated RnBeads output) 
2. [Raw beta value and detection p-value extraction](scripts/preprocessingSEst)
    + [sex prediction using sEst](scripts/sexPrediction.R)
    + [visualisation of the predicted sample sexes](scripts/sexPredictionVisualisation)
    
## Principal component analysis

1. [Calculation of principal components](scripts/)
    + [association of PCs with sample characteristics & heatmap visualisation](scripts/)
  
## Cell composition deconvolution

1. [Preprocessing for cell composition prediction (MethylumiNoob method)](scripts/preprocessingMinfiNoob)
    + [placental lineage cell composition prediction](scripts/)
    + [visualisation of the results](scripts/)

## Lineage-specific methylation characterisation
[Additional scripts](otherScripts/) to look at methylation features that are specific to placental / fetal lineages are provided. Supporting data required for these analysis can be found [here](supportingData).

## References
Reference datasets were obtained from the following publications:





