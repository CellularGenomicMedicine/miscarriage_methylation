# miscarriage_methylation
Methylation analysis of EM and CV samples from first trimester pregnancy loss to conduct cell composition deconvolution and to determine whether the sampled tissue types are distinct and whether they correspond to the expected cell sample types. The methylome analyses were carried out on data generated using the EPIC array.

## Pre-processing & sample QC
1. [SWAN normalisation and sample / site QC](scripts/preprocessingSWAN.R)
    (the SNP heatmap can be extracted from the generated RnBeads output) 
2. [Raw beta value and detection p-value extraction](scripts/preprocessingsEst.R)
    + [sex prediction using sEst](scripts/sexPrediction.R)
    + [visualisation of the predicted sample sexes](scripts/sexPredictionVisualisation.R)
    
## Principal component analysis

1. [Calculation of principal components](scripts/principalComponentAnalysis.R)
    + association of PCs with sample characteristics & heatmap visualisation
  
## Cell composition deconvolution

1. [Preprocessing for cell composition prediction (MethylumiNoob method)](scripts/preprocessingMinfiNoob.R)
    + [placental lineage cell composition prediction](scripts/cellDeconvolution.R)
    + [visualisation of the results](scripts/cellCompositionVisualisation.R)

## Lineage-specific methylation characterisation
[Additional scripts](otherScripts/) to look at methylation features that are specific to placental / fetal lineages are provided. Supporting data required for these analysis can be found [here](supportingData). The supporting data is as follows:
1. EPIC sites required for cell deconvolution of first trimester placental samples (Yuan et al.)
2. EPIC sites within repetitive Line1 and Alu regions
3. EPIC sites within placenta-specific and generic imprinted regions
4. EPIC sites within placenta-specific partially methylated domains

## References
1. Yuan, V., Hui, D., Yin, Y. et al. Cell-specific characterization of the placental methylome. BMC Genomics 22, 6 (2021). https://doi.org/10.1186/s12864-020-07186-6
2. Schroeder DI, Blair JD, Lott P, Yu HOK, Hong D, Crary F, et al. The human placenta methylome. Proc Natl Acad Sci. 2013;110(15):6037–42





