---
title: "PCRedux: machine learning helper tool for sigmoid curves"
tags:
  - "R"
  - "PCR"
  - "quantitative Polymerase Chain Reaction"
  - "machine learning"
  - "sigmoid"
authors:
  - name: "Stefan Rödiger"
    orcid: 0000-0002-1441-6512
    affiliation: 1
  - name: "Michał Burdukiewicz"
    orcid: 0000-0001-8926-582X
    affiliation: 2
  - name: "Roy-Arne Senkel"
    orcid: 0000-0001-6139-9471
    affiliation: 1
  - name: "Andrej-Nikolai Spiess"
    orcid: 0000-0002-9630-4724
    affiliation: 3
  - name: "Thomas Fischer"
    orcid: 0000-0001-6235-2261
    affiliation: 1
affiliations:
  - name: "Brandenburg University of Technology Cottbus-Senftenberg, Germany"
    index: 1
  - name: "Uniwersytet Wrocławski: Wrocław, Poland"
    index: 2
  - name: "University Medical Center Hamburg-Eppendorf, Hamburg, Germany"
    index: 3
date: "31 October 2017"
bibliography: literature.bib
output:
  html_document:
    keep_md: yes
---

# Summary

Data having a sigmoid ('S'-shaped) curves are found in many experiments. An 
example are amplification curve data from quantitative Polyerase Chain Reactions 
(qPCR). The qPCR is an indispensable method in many disciplines such as human 
diagnostics and forensics [@martins_dna_2015]. There are numerous software 
packages for the analysis of qPCR data were published [@pabinger_survey_2014, 
@rodiger_r_2015]. These software packages provide complex pipelines and a rich 
sets of criteria to process qPCR data adequately. This includes the 
pre-processing of raw data, fitting of non-linear models on raw data, the 
calculation of the cycle of quantification the calculation of amplification 
efficiencies, the relative gene expression analysis, normalization procedures, 
and data management [@rodiger_chippcr:_2015, @spiess_impact_2015, 
@rodiger_enabling_2017, @mallona_chainy:_nodate].

The `PCRedux` package contains functions and data sets for machine learning and 
statistical analysis with a focus on amplification curve data. The data sets 
originate form qPCRs, which were created and rated by human experts. 
Amplification curves have characteristics, which can be used for the 
classification. Analysis on the curve data such as changepoint analysis, 
regression analysis, noise analysis, autocorrelation analysis and model fitting 
are applied to the qPCR data and yield more than 30 feautures. This can be used 
for the creation of models that predict a class (e.g., positive, ambiguous, 
negative qPCR reaction) from input features (slope, background level, 
changepoints) based on implementations by others (e.g., @erdman_bcp:_2007, 
@Ritz2008, @Febrero_Bande_2012, @james_ecp:_2013) and us [@roediger_RJ_2013, 
@rodiger_chippcr:_2015]. The features from amplification curves can be extract 
by the `pcrfit_helper` function. `pcrfit_parallel` is a wrapper for numerous 
functions, which are performed in parallel. Since this function provides a 
complex set of features we included the visualization function 
`visdat_pcrfit_parallel`, which is based on the `vis_dat` function by 
@Tierney2017. In addition, the following function were developed:

- `autocorrelation_test` performs an autocorrelation analysis on qPCR data,
- `decision_modus` finds the most frequent rating by a human,
- `earlyreg` performs a regression analysis in background region,
- `head2tailratio` compares the ratio of the head and tail,
- `hookregNL` and `hookreg` attempt to detect a hook effect in the amplification curve,
- `mblrr`, a function to perform a qunantile-filter based local robust regression,
- `performance`, performance analysis (e.g., sensitivity, specificity, Cohen's kappa) for binary classification, and
- `qPCR2fdata`, a helper function to convert amplification curve data to the *fdata* format.



`PCRedux` provides tools for features extraction from sigmoid data and machine 
learning. The `PCRedux` package an add-on package for the open source 
statistical computing language and environment *R* [@R_language].

![](fig1.png)<!-- -->

# Acknowlegements
This work was funded by the Federal Ministry of Education and Research
(BMBF) InnoProfile-Transfer-Projekt 03IPT611X and in part by 'digilog: Digitale
und analoge Begleiter für eine alternde Bevölkerung' (Gesundheitscampus
Brandenburg, Brandenburg Ministry for Science, Research and Culture).

Corresponding Author: stefan.roediger@b-tu.de

\newpage

# References
