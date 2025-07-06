# depression_transcriptome

This code corresponds to:

Systematic review and mega-analysis of the peripheral blood transcriptome in depression implicates dysregulation of lymphoid cells and histones

Authors: Chaitanya Erady, Richard A.I. Bethlehem, Ed Bullmore, Mary-Ellen Lynall

## Files:

### Differential gene and transcript mega-analyses
* __xxx_data_processing.Rmd__ Data processing instructions for individual datasets
* __xxx_DGE.Rmd__ DGE analysis of individual datasets, Figure S4
* __template_for_DTE.Rmd__ DTE analysis of individual datasets
* __meta_analysis_with_bias_and_inflation_correction.Rmd__ Conducts bias and inflation corrected meta-analysis: DGE cell corrected, DGE not cell corrected, BMI sensitivity analysis, Platelet count sensitivity analysis, DTE cell corrected meta-analysis
* __gene_venn.Rmd__ Creates Venn diagrams, Figure 1A
* __qq_plots.r__ Creates qq plots for meta-analysed results, Figure S5
* __LOO_analyses.Rmd__ Conducts leave-one-out DGE analysis, Data for Table S8
* __bacon_weightedz.r__ Conducts bias and inflation corrected weighted Z meta-analysis
* __bacon_weightedz_withNAs.r__ Conducts bias and inflation corrected weighted Z meta-analysis when data has NA values
* __wZ.r__ Code to conduct weighted Z meta-analysis

### Differential transcript usage
* __dtu_function.r__ Functions needed for DTU analysis
* __template_for_DTU.Rmd__ DTU analysis of individual datasets

### Comparisons
* __compare_DGE_with_proteome.Rmd__ Compares mega-analytic DGE results to differential protein expression from Daskalakis 2024 _Science 384:eadh3707_
* __RRHO_plots.Rmd__ Compares DGE to whole blood TWAS results from Meng X et al. 2024 Nat Gen, generates RRHO plots, Figures 1 B,C, and S6 B, C 



### Core cellular processes mega-analysis
* __gene_group_analyses.Rmd__ Creates data for gene groups correspondign to core cellular processes, permutation test to assess significance, meta-analysis of permutation results
* __gene_group_perm_plot.Rmd__ Creates Figures 2C, S8


### Gene/transcript enrichment and cell origin analyses
* __GSEA_plots.Rmd__ Conducts enrichment analysis using Reactome and MSigDB, create leading edge plots, Figures 2A, S7
* __gsea_leading_edge.r__ Code to create leading edge plots
* __dge_sc_enrichment.Rmd__ Performs and plots LR Cell analyses, Figures 2D, S9
* __plot_trynka_enrichment.Rmd__ Makes T cell activation enrichment plot, Figure 2B



### WGCNA meta-analysis
* __consensus_clustering.Rmd__ Conducts consensus WGCNA for CNT+MDD samples, Figure 3A,B
* __consensus_clustering_control_samples_only.Rmd__ Conducts consensus WGCNA for CNT samples only
* __WGCNA_module_exploration.Rmd__ Codes for analysis of module enrichment, module preservation scores, module-trait association, module-disorder association, identify hub genes, module membership of DGE genes, inter-dataset correlations, module-smoking association, leave-one-out analysis: without dbGaP, calculate average weighted cohensd, Figures 6A, S10, S11, Data for Tables S4, S5, S6 and S7
* __wgcna_cohensd.Rmd__ Makes WGCNA meta-analysis plots, Figures 3C,D


## Processing additional MDD case-control datasets
* To meta-analyse additional MDD case-control datasets, see steps listed for BIODEP/Le/Mostafavi data processing if using RNA-Seq data, else see steps listed for dbGaP/HiTDiP data processing if using microarray data.
* Use meta_analysis_with_bias_and_inflation_correction.Rmd to conduct meta-analysis of processed datasets.
* Use template_for_DTE.Rmd and template_for_DTU.Rmd, respectively, to conduct DTE and DTU analysis.
* Gene-, transcript-level summary statistics and harmonised processed individual-level count matrices along with estimated cell counts and metadata are provided (where permitted) at the Zenodo repository https://doi.org/10.5281/zenodo.15290507. 

