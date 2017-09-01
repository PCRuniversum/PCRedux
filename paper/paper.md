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
affiliations:
  - name: "Brandenburg University of Technology Cottbus-Senftenberg, Germany"
    index: 1
  - name: "Uniwersytet Wrocławski: Wrocław, Poland"
    index: 2
date: "19 Septemer 2017"
bibliography: paper.bib
output:
  html_document:
    keep_md: yes
---

# Summary

Data having a characteristic sigmoid ('S'-shaped) curves are commonly found in 
experimental data. An example are amplification curve data from quantitative 
Polyerase Chain Reaction (qPCR) experiments. There are numerous software 
packages for R for the analysis of quantitative PCR (qPCR) data 
[@pabinger_survey_2014, @rodiger_r_2015]. so far, there was no R package, which 
prepares amplification curve data or extracts curve features for machine 
learning.

## General workflow

The `PCRedux` package contains function and data sets for machine learning on 
quantitative PCR data in R. They are intended for the creation of models that 
predict a class (positive, ambigous, negative qPCR reaction) from input features 
(slope, hook effect, background level, changepoints).

### High-level functionality and purpose of the software

Amplification curve characteristics such a the hook effect are challenging 
during the amplification curve classification. Therefore, the unique hookreg 
function was introduced to determine if an amplification curve exhibits this 
characteristics.

### Data sets

The data sets encompass amplification curve data for training purposes.

### Use cases

- A summary describing the high-level functionality and purpose of the software
for a diverse, non-specialist audience
- A clear statement of need that illustrates the purpose of the software
- A list of key references including a link to the software archive

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

![](hookreg.png)<!-- -->

# Acknowlegements
This work was funded by the Federal Ministry of Education and Research
(BMBF) InnoProfile-Transfer-Projekt 03IPT611X and in part by 'digilog: Digitale
und analoge Begleiter für eine alternde Bevölkerung' (Gesundheitscampus
Brandenburg, Brandenburg Ministry for Science, Research and Culture).

\newpage

# References
