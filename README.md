# EAGLE_LUAD
Genetic and Epigenetic Intratumor Heterogeneity Impacts Prognosis of Lung Adenocarcinoma


This resource provides the code developed in the study of Hua et al. "Genetic and Epigenetic Intratumor Heterogeneity Impacts Prognosis of Lung Adenocarcinoma". It reproduces the key results of the survival analysis using variance-weighted Cox proportional-hazards model, and can be applied to other genomic studies with multi-sampling design to explore the association of intratumor heterogeneity and survival.

Requirements:
R (tested in R version 3.6.0)
R libraries: survival, survey

Data:
The data for survival analysis in this study is provided in the package (/data/data4demo.RData).

Quick start:
To reproduce the result reported in Hua et al., unzip the package. In R go to the /code directory and run source(code4Demo.R). (Estimated running time < 1min)

General notes:
The code provided in code4Demo.R reproduces the survival analysis using variance-weighted Cox proportional-hazards model. It also generates the study figures and table in the directory of /figures and /table.

Citation:
Hua et al. "Genetic and Epigenetic Intratumor Heterogeneity Impacts Prognosis of Lung Adenocarcinoma".
