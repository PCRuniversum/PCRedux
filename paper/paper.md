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
  - name: "University  Medical  Center  Hamburg-Eppendorf,  Hamburg, Germany"
    index: 3
date: "19 Septemer 2017"
bibliography: paper.bib
output:
  html_document:
    keep_md: yes
---

# Summary

Data having a characteristic sigmoid ('S'-shaped) curves are found in many 
experiments. An example are amplification curve data from quantitative Polyerase 
Chain Reaction (qPCR) experiments. There are numerous software packages for for 
the analysis of qPCR data. In particular, software packages for the statistical 
computing language R were published [@pabinger_survey_2014, @rodiger_r_2015]. so 
far, there was no R package, which prepares amplification curve data or extracts 
curve features for machine learning. `PCRedux` provides tools for importing and 
working with sigmoid data.

The `PCRedux` package contains functions and data sets for machine learning and 
statistical analysis with a focus on sigmoid (amplification curve) data. In 
detail, package contains data sets of qPCR, which were created and rated by 
human experts. Amplification curves have characteristics which can be used for 
the classification. The features from amplification curves can be extract by the 
`pcrfit_parallel` function. `pcrfit_parallel` performs in parallel multiple 
analysis on the curve data such as changepoint analysis, regression analysis, 
noise analysis, autocorrelation analysis and model fitting to qPCR data. They 
can be used for the creation of models that predict a class (e.g., positive, 
ambiguous, negative qPCR reaction) from input features (slope, background level, 
changepoints) based on implementations by others (e.g., @erdman_bcp:_2007, 
@ritz_qpcr:_2008, @Febrero_Bande_2012, @james_ecp:_2013) and us 
[@rodiger_surface_2013, @rodiger_chippcr:_2015]. 

![](fig1.png)<!-- -->

# Acknowlegements
This work was funded by the Federal Ministry of Education and Research
(BMBF) InnoProfile-Transfer-Projekt 03IPT611X and in part by 'digilog: Digitale
und analoge Begleiter für eine alternde Bevölkerung' (Gesundheitscampus
Brandenburg, Brandenburg Ministry for Science, Research and Culture).

Corresponding Author: stefan.roediger@b-tu.de

\newpage

# References
