# -rupp_tomlinson-jclininvest-2023
# Rupp et al. 2023

All the data and scripts necessary to generate the analysis and figures presented in Rupp et al. mansucript (see citation below).

## Files
```
+- data
    +- GSE162603
        +- GSE162603_counts.csv.gz
        +- GSE162603_metadata.csv.gz
    +- GSE172203
        +- GSE172203_barcodes.tsv.gz
        +- GSE172203_features.tsv.gz
        +- GSE172203_matrix.mtx.gz
        +- GSE172203_metadata.csv.gz
    +- mbh
        +- NovaA-219
        +- NovaA-225
        +- NovaA-230
    +- sun1
        +- NovaA-296
    +- 
+- scripts
    +- conservation.R
    +- figure_functions.R
    +- figures.R
    +- GSE162603.R
    +- scrublet.py
    +- snrna-seq_functions.R
    +- sun1.R
+- environment.yml
+- README.md
```

## Pipeline
This pipeline is best run with [anaconda](https://www.anaconda.com/products/individual) to replicate the original environment. This is optional, assuming all available packages and tools are installed and the appropriate version. See `environment.yml` for the package versions used in this analysis.
### Set up `conda` environment
1. Create the `conda` environment: `conda env create -f environment.yml`
2. Activate the `conda` environment: `conda activate Glp1r`
### Analysis
1. Analyze published Lepr TRAP-seq data: `Rscript scripts/GSE162603.R`
2. Analyze the Lepr-Sun1 snRNA-seq data: `Rscript scripts/sun1.R`
3. Analyze the combined mouse/rat/macaque hypothalamus snRNA-seq data: `Rscript scripts/conservation.R`
4. Generate the figure panels: `Rscript scripts/figures.R`

## Citation
Rupp AC, Tomlinson AJ, Affinati AH, Yacawych WT, Duensing AM, True C, Lindsley SR, Kirigiti MA, MacKenzie AJ, Polex-Wolf J, Li C, Knudsen LB, Seeley RJ, Olson DP, Kievit P, Myers MG Jr. Suppression of food intake by Glp1r/Lepr-coexpressing neurons prevents obesity in mouse models. J Clin Invest. 2023 Aug 15:e157515. doi: 10.1172/JCI157515. Epub ahead of print. PMID: 37581939.

