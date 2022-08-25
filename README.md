# Master Thesis

The order in which to run the needed code is shown below

## Alignment

1. `alignment/indexation.R`: Build index for transcriptomic analysis

2. `alignment/alignment_ncbi.R`: Align the reads from FastQ files

3. `alignment/counts_ncbi.R`: Generate count matrix from the .bam files

## Differential Expression Analysis

4. `analysis_counts/analysis_counts.R`: DEG analysis with DESeq2, limma and edgeR.

5. `analysis_counts/gene_intersection.R`: Get gene intersection between the DEGs obtained with the three methods.

## Functional Analysis

6. `ClusterProfiler/Cluster_profiler.R`: Apply ORA and GSEA for functional analysis based on the Msigdb gene sets.

## Tumor Characterization

7. `Characterization/GSVA_mouse.R`: Use GSVA to obtain an enrichment score profile for the selected categories.

8. `Characterization/IL10_differences.R`: Extract the gene with the greatest logFC change from the Interleukin 10 pathway.

9. `Characterization/convert_mouse_to_human.R`: Convert mice gene symbols in the categories from mice to human.

10. `Characterization/human_GSVA.R`: Use GSVA and the translated categories to obtain a profile for the human samples.

11. `Characterization/Tumor_classification.ipynb`: Use SVM to classify the human samples as pBIC-like and pBIC10-like.

12. `Characterization/categorization_analysis.R`: Analize the categorized ssamples to find characteristic genes that differentiate between pBIC-like and pBIC10-like

## Other Files

The following files are not required to proceed with the analysis

- `ClusterProfiler/exploratory_plots.R`

- `alignment/multiqc_plot.R`

- `analysis_counts/visualization.R`
