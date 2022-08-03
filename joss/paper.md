---
title: 'PCRedux: A Quantitative PCR Machine Learning Toolkit'
tags:
  - R
  - qPCR
  - bioinformatics
  - machine learning
  - reproducible research
date: "12 May 2022"
affiliations:
  - name: Autonomous University of Barcelona, Bellaterra, Spain
    index: 1
  - name: Medical University of Białystok, Białystok, Poland
    index: 2
  - name: Soilytix GmbH, Hamburg, Germany
    index: 3
  - name: Warsaw University of Technology, Warsaw, Poland
    index: 4
  - name: Pirogov Russian National Research Medical University, Moscow, Russia
    index: 5
  - name: BTU Cottbus–Senftenberg, Faculty of Health Brandenburg, Senftenberg, Germany
    index: 6
  - name: BTU Cottbus–Senftenberg, Faculty Environment and Natural Sciences, Senftenberg, Germany
    index: 7
authors:
  - name: Michał Burdukiewicz^[Co-first author]
    orcid: 0000-0001-8926-582X
    affiliation: 1, 2
  - name: Andrej-Nikolai Spiess^[Co-first author]
    orcid: 0000-0002-9630-4724
    affiliation: 3
  - name: Dominik Rafacz
    orcid: 0000-0003-0925-1909
    affiliation: 4
  - name: Konstantin Blagodatskikh
    orcid: 0000-0002-8732-0300
    affiliation: 5
  - name: Stefan Rödiger^[Corresponding author]
    orcid: 0000-0002-1441-6512
    affiliation: 6, 7
bibliography: paper.bib
---

# Summary

qPCR (quantitative polymerase chain reaction) is indispensable in research,
diagnostics and forensics, because it provides quantitative information
about the amount of DNA in a sample [@pabinger_survey_2014]. The interpretation
of amplification curves (ACs) is often difficult if the curve does not follow a typical sigmoidal trajectory.

`PCRedux` is an `R` package [@stats] for feature extraction and classification
in the realm of explainable machine learning, which uses statistical functions to compute
90 boolean and numerical descriptors from ACs. It can also be used to determine *Cq* values and amplification
efficiencies (*E*) for high-throughput analysis.

Given the lack of class-labeled qPCR data sets, `PCRedux` includes functions for
aggregation, management and dissemination of qPCR
datasets that can, but must not necessarily be, trichotomously classified into _negative_,
_positive_ and _ambiguous_ curves.

# Statement of need

qPCR is a widely used laboratory method for the precise detection and quantification of pathogens and gene expression.
The latter has contributed significantly to the understanding of physiological and pathological processes in pharmacology, medicine and forensics.
[@pabinger_survey_2014; @kok_small_2018] Although available software packages provide workflows and criteria for processing qPCR
data (pre-processing of raw data, fitting of non-linear models, calculation of a threshold- or second derivative-based *Cq* or *E*, relative gene expression 
analysis, normalization procedures and data management), they lack functionality for machine learning.
[@pabinger_survey_2014; @ruijter_evaluation_2013; @ruijter_efficiency_2021;
@ramakers_assumption-free_2003].

qPCR curves must meet quality criteria for analysis and are often categorized 
by the user according to rather subjective criteria (_e.g._,
sigmoidal shape, slope, noise, presence of a "hook effect")
[@burdukiewicz_algorithms_2018; @spiess_system-specific_2016;
@spiess_impact_2015; @hanschmann_looptag_2021_2]. While positive qPCR reactions 
usually exhibit a sigmoidal shape, negative ACs display a rather flat and linear trajectory (\autoref{fig:fig_1}).

![Analysis of ACs using the `PCRedux` package. A) ACs exhibit a high diversity in their appearance. The left plot (positive)
shows ACs of which almost all are sigmoidal. The
signal amplitude ranges from - 70 to 4000 relative fluorescence units (RFU).
Of importance are those ACs that go slightly into the negative
range of 10 - 20 cycles. Top right (negative) are negative
ACs (signal amplitude between - 700 and 300) with noise (cycle
20 - 40). The bottom right plot (ambiguous) shows ACs that do
not possess typical nonlinear slopes (non-sigmoid). These cannot be clearly
classified as positive or negative. B) From A), the values of the three descriptors cpD2 (maximum of the second
derivative, equals *Cq*), amptester_polygon (area under the curve) and cpDdiff (absolute
cycle distance between the maximum of the second and first derivative) were
calculated and plotted for the three classes. Data
from htPCR dataset [@ritz_qpcr:_2008].\label{fig:fig_1}](fig_1.png)

So how can ACs be objectively and
reproducibly assessed and automatically interpreted (_e.g._, as
positive/negative/ambiguous or low/high quality)? For high-throughput
experiments, manual evaluation is not feasible because of mental exhaustion errors or non-reproducibility from arbitrary thresholds or subjective assessments. 
While internal laboratory guidelines seem to partially remedy this, they are usually not standardized with other labs.
[@bustin_why_2010; @taylor_ultimate_2019; @kim_experimenting_2018].

Automatically extracted features from ACs (_e.g._, *Cq* and *E*, slopes, change points, features of local curve segments) could provide a
solid feature basis for classification by machine learning. Yet to date, no open-source
software applies classical biostatistical methods for explainable machine
learning on ACs. Furthermore, there are no class-labelled
datasets that can be used in this context.

`PCRedux` is the first open-source software that can extract 90 mathematical
descriptors (features) from raw ACs. The features are
numerically or analytically derived, quantifiable, informative properties of scaled ACs in scalar units.

# Software engineering

`PCRedux` (v.1.1-2, [MIT license](https://mit-license.org/)) is an `R` package
(S3 class system). `R` was chosen because it provides comprehensive tools for reproducible
statistical and bioinformatics analyses [@gentleman_bioconductor:_2004;
@gentleman_statistical_2007; @rodiger_r_2015; @liu_r_2014;
@leeper_archiving_2014]. Unit tests using the `testthat` package [@testthat]
were used for software quality control of `PCRedux`.

## Functions

Conceptually, we divide ACs into regions of interest (ROI) for feature
calculation (vignette @PCRedux Figure 5). Typical for qPCR, baseline, exponential/linear and plateau phases are
located at the left, middle and right tail region of the curve, respectively (vignette @PCRedux Figure 5).

`PCRedux`'s algorithms, published by others and ourselves (`qpcR`
[@ritz_qpcr:_2008], `MBmca` [@rodiger_surface_2013], `chipPCR`
[@rodiger_chippcr:_2015]), were adapted for qPCR analysis. The `PCRedux`
dependencies include packages for preprocessing (`chipPCR`,
`MBmca`), fitting of non-linear models and calculation of
*Cq* and *E* (`qpcR`). `pcrfit_single()` is the workhorse function for single ACs that generates a *data.frame* with 90 descriptors. 
All output values are of type `numeric`, even if boolean. Among others, we included autocorrelation analyses, (Bayesian) change-point analyses (`bcp`
[@erdman_bcp:_2007], `ecp` [@james_ecp_2015]), area determinations (`pracma`
[@pracma]), regression (multi-parametric non-linear & robust local & regression
models with segmented relationships: (`robustbase`
[@todorov_object-oriented_2009], `stats` [@stats], `segmented`
[@muggeo_interval_2017]) and hook effect detection (`PCRedux`
[@burdukiewicz_algorithms_2018]).

`encu()` is an extension of `pcrfit_single()`, where meta
information such as detection chemistry and platform is included, and is suitable
for large data sets. Both functions are error-proof and utilize, among others, the following
descriptor-generating functions:

* `earlyreg()`, calculates features by regression analysis in the background region,
* `head2tailratio()`, compares the ratio of head and tail,
* `hookregNL()` and `hookreg()`, try to detect a hook effect [@burdukiewicz_algorithms_2018],
* `mblrr()`, performs a local robust regression analysis,
* `winklR()`, calculates the angle based on the first, and the second derivative and
* `autocorrelation_test()`, tests for autocorrelation.

Auxiliary preprocessing and analysis functions of the package are:

* `armor()`, catches errors and creates the output,
* `decision_mode()`, calculates the frequency of classes in a dataset,
* `qPCR2fdata()`, converts AC data to the *fdata* format for Hausdorff distance analysis [@fda.usc] and
* `performeR()`, performs power analyses (e.g., sensitivity, specificity, Cohen's $\kappa$) for binary classification.

Application examples in the context of machine learning can be found in the @PCRedux vignette.

## Graphical User Interface:

`run_PCRedux()` invokes a graphical user interface (\autoref{fig:fig_2}) based
on the `Shiny` technology [@shiny], providing features as a downstream accessible table.

![Graphical user interface for the analyses of qPCR data. A) The `run_PCRedux()`
GUI for analysis and tabular display can use browsers or R environments that
support `ECMA Script` and `HTML`. In this example, the GUI was used in `RKWard`
(v.0.7.2, Linux, Kubuntu 21.10, [@rodiger_rkward:_2012]). B) Optionally, information about the current state of errors can be obtained via the R console.\label{fig:fig_2}](fig_2.png)

## Datasets and Data Labeling

`PCRedux` contains class-labeled ACs (n = 14360; label: negative,
positive, ambiguous) from various qPCR instruments and detection methods, as
determined by the majority vote of four experienced researchers (@PCRedux). 
Class labels were derived from the `humanrater2()` function, which uses `tReem()` for shape similarity calculation [@fda.usc].

# Conclusion

Manual classification of ACs is time-consuming and error-prone,
especially with large data sets where a significant proportion of curves deviate from sigmoidal shape,
and where the results are influenced by subjective perception. An automated system for analyzing qPCR curves
offers objectification and generalization of the decision-making process.

Here, training neural networks poses a viable option. The question is what the
resulting network considers relevant, especially in a diagnostic scenario. 
To avoid these 'black box' situations, `PCRedux` enables a fast and computer-assisted classification of ACs based on 90
statistically and analytically founded descriptors, aiming to improve the quality and reproducibility of qPCR data analysis.

# Acknowledgments

Grateful thanks belong to the `R` community.

# Funding

None

# References
