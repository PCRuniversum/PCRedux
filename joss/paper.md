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
- name: Medical University of Białystok, Białystok, Poland
  index: 1
- name: Soilytix GmbH, Hamburg, Germany
  index: 2
- name: Warsaw University of Technology, Warsaw, Poland
  index: 3
- name: Pirogov Russian National Research Medical University, Moscow, Russia
  index: 4
- name: BTU Cottbus–Senftenberg, Faculty of Health Brandenburg, Senftenberg, Germany
  index: 5
- name: BTU Cottbus–Senftenberg, Faculty Environment and Natural Sciences, Senftenberg,
    Germany
  index: 6
authors:
- name: Michał Burdukiewicz^[Co-first author]
  orcid: 0000-0001-8926-582X
  affiliation: 1
- name: Andrej-Nikolai Spiess^[Co-first author]
  orcid: 0000-0002-9630-4724
  affiliation: 2
- name: Dominik Rafacz
  orcid: 0000-0003-0925-1909
  affiliation: 3
- name: Konstantin Blagodatskikh
  orcid: 0000-0002-8732-0300
  affiliation: 4
- name: Stefan Rödiger^[Corresponding author]
  orcid: 0000-0002-1441-6512
  affiliation: 5, 6
bibliography: paper.bib
---

# Summary

Numerous scientific data with sigmoidal ('S'-shaped) curves exist, with one of
the more important use cases being amplification curve data from quantitative
polymerase chain reactions (qPCR). qPCR has developed into an indispensable
method in research, human diagnostics and forensics, because quantitative
information can be determined from an amplification curve that provides
information about the amount of DNA in a sample. An important descriptor here is
the 'cycle of quantification' (*Cq*) which can be used for gene expression
analysis or viral load determination. Software packages are available that
provide workflows and criteria for processing qPCR data, including preprocessing
of raw data, fitting of non-linear models to raw amplification curve data,
calculation of threshold or second-derivative based *Cq* values, calculation of
amplification efficiency, relative gene expression analysis, normalization
procedures and data management [@pabinger_survey_2014]. At the same time, the
interpretation of amplification curves can often be challenging if the curve
trajectory does not follow a typical amplification curve, for instance how to
objectively assess the quality of an amplification curve or how to automatically
interpret amplification curves as positive, negative, or ambiguous.

To date, there is no open-source software package that provides classical
biostatistical methods-based ("explainable") machine learning on amplification
curves. Furthermore, no class-labeled datasets exists that can be used in this
context.

The `PCRedux` package is an add-on package ([MIT
license](https://mit-license.org/)) for the open-source statistical programming
language and environment `R` [@R_Core_Team]. This package contains functions
that determine features from amplification curves using classical statistical
and analytical procedures used in the classification of amplification curves.
These 90 procedures are Boolean and numeric descriptors that can be employed for
downstream machine learning and classification. In addition, the `PCRedux`
package can be used for determination of *Cq*-values and amplification
efficiencies (*E*) in high-throughput analysis. Finally, `PCRedux` provides an
extensive cohort of labeled amplification curve datasets obtained from various
qPCR instruments and detection methods, where the amplification curves were
trichotomously classified by several qPCR users into _negative_, _positive_ and
_ambiguous_.

In summary, the `PCRedux` package is a versatile tool for feature extraction and
classification of amplification curves based on explainable machine learning.

# Statement of need

Quantitative real-time PCR (qPCR) is a widely used laboratory method for the
precise quantification of nucleic acids, *e.g.*, in the detection and
quantification of pathogens such as viruses and bacteria. Another application of
qPCR is gene expression analysis, which has contributed significantly to our
understanding of human physiological and pathological processes in pharmacology,
medicine and forensics. [@pabinger_survey_2014; @kok_small_2018]. Typically, the
workflow of qPCR data analysis includes the determination of the cycle of
quantification (*Cq*) and the amplification efficiency (*E*)
[@ruijter_evaluation_2013; @ruijter_efficiency_2021;
@ramakers_assumption-free_2003].

In order to perform this analysis, the amplification curves must fulfill certain
(quality) criteria. Amplification curves can be categorized by the user on the
basis of somewhat subjective criteria, such as sigmoidal shape, steepness,
noisiness and overall "quality". Positive qPCR reactions usually exhibit a
sigmoidal shape, while negative amplification curves will have a more flat and
linear appearance [@ritz_qpcr:_2008]. The quality of the amplification reaction
can be identified by - among other things - whether they are noisy, include a
"hook effect", or show other disturbances [@burdukiewicz_algorithms_2018;
@spiess_system-specific_2016; @spiess_impact_2015; @hanschmann_looptag_2021_2].
Figure \autoref{fig:fig_1} A shows an experiment with several amplification
curves that exhibit different curvatures.

![Analysis of amplification curves using the `PCRedux` package. A) Amplification
curves can have a high diversity in their appearance. The left plot (positive)
shows several amplification curves, almost all of which are sigmoidal. The
signal amplitude ranges from - 70 to 4000 relative fluorescence units (RFU).
Noticeable are the amplification curves, which go slightly into the negative
range in the range of 10 - 20 cycles. Top right (negative) are negative
amplification curves (signal amplitude between - 700 and 300) with noise (cycle
20 - 40). The bottom right plot (ambiguous) shows amplification curves that do
not have typical nonlinear slopes (non-sigmoid). They cannot be clearly
classified as positive or negative. B) From the data from A, the values (Value,
[A.U.] (arbitrary unit)) of three descriptors cpD2 (maximum of the second
derivative), amptester_polygon (area under the curve) and cpDdiff (absolute
cycle distance between the maximum of the second and first derivative) were
calculated and plotted for the classes positive negative and ambiguous. Data
from htPCR dataset [@ritz_qpcr:_2008].\label{fig:fig_1}](fig_1.png)

However, the question arises on how the quality of an amplification curve can be
objectively assessed and on how amplification curves can be automatically
interpreted as positive, negative or ambiguous in a reproducible manner, i.e.
excluding the human factor. Especially in high-throughput experiments, manual
analysis of all amplification curves is not readily feasible. First, it is
highly tedious and becomes error-prone due to exhaustion of concentration
ability. Second, manual analyses are often not reproducible as scientists in
different laboratories use different manual thresholds or incorporate a
subjective evaluation of the amplification curve quality. Here, in-house
laboratory guidelines appear to facilitate the evaluation of amplification
curves, but these do not transfer to other laboratories [@bustin_why_2010;
@taylor_ultimate_2019; @kim_experimenting_2018]. Consequently, and regardless of
sample throughput, these amplification curves must be classified quickly
according to objective criteria to determine whether they pertain to one of the
three (or more) classes mentioned above.

To solve this problem, automatically extracted features from amplification
curves could provide a solid ground for the classification of amplification
curves by machine learning. Numerical descriptors are features such as *Cq* and
"*E*" values, slopes, change points, features of local curve segments and others
that differ between amplification curves. Training more or less "black box"
neural networks would pose an option, but here the question of what the model
considers as relevant would arise, especially in diagnostic applications. In
order to avoid these situations, methods are available in which each descriptor
can be mathematically explained. So far, and to the best of our knowledge, there
is no open-source software that provides an automatic and reproducible
classification of amplification curves based on statistical and mathematically
determined machine learning features.

`PCRedux` is the first open-source software to automatically determine 90
mathematical descriptors (features) extracted from raw amplification curves that
can be subjected to downstream machine learning algorithms, where a feature is a
numerically or analytically derived quantifiable informative property of an
amplification curve in a scalar unit. Consequently, the software can be used as
a feature-delivering step before an already existing machine learning-based
classifications pipeline.

Scientific work always delivers and depends on data. In particular, open data
are the cornerstones of science. Given the lack of datasets with classified
amplification curves, the development of the `PCRedux` package also addressed
the aggregation, management and distribution of classified qPCR datasets. This
involved the aggregation of more than 14,000 amplification curves obtained from
various thermocyclers and detection systems. In addition, it includes a function
for blinded random classification of amplification curves.

`PCRedux` is designed to improve the quality and reproducibility of qPCR data
analysis by systematically validating input data. It also supports the
development of control mechanisms that can be used for other automated
algorithms (e.g., *Cq* value determination) and provides open data. We hope that
our results will accelerate the automation of qPCR data analysis. In the optimal
case, `PCRedux` can reduce the workload in analyzing large amounts of data and
provide reproducible and more objective data analysis.

# Software engineering

`PCRedux` (v.\~1.1-2) is an `R` package (S3 class system) that calculates
features from amplification curve data of quantitative polymerase chain
reactions (qPCR) for machine learning purposes. `R` was chosen because it
provides comprehensive tools for reproducible statistical and bioinformatics
analyses [@gentleman_bioconductor:_2004; @gentleman_statistical_2007;
@rodiger_r_2015; @liu_r_2014; @leeper_archiving_2014]. Unit tests using the
`testthat` package [@testthat] were used for software quality control of
`PCRedux`.

## Functions

An amplification curve is divided into regions of interest (ROI) for feature
calculation (vignette @PCRedux Figure 5). For example, the baseline phase is
located in the left tail region, while the plateau phase is located in the right
tail region. The exponential and linear phase is located between these two ROIs
(vignette @PCRedux Figure 5).

We integrated into `PCRedux` algorithms adapted for qPCR analysis published by
others and ourselves (`qpcR` [@ritz_qpcr:_2008], `MBmca`
[@rodiger_surface_2013], `chipPCR` [@rodiger_chippcr:_2015]). `pcrfit_single()`
is the main workhorse function for descriptor calculation of single
amplification curves. The `PCRedux` dependencies thus include packages that
contain workflows and criteria for processing qPCR data. These include
preprocessing of raw data (`chipPCR`, `MBmca`), fitting of non-linear models to
raw data, calculation of quantification cycles and amplification efficiencies
(`qpcR`). The function generates a *data.frame* with 90 descriptors from an
amplification curve trajectory. All output values are of type `numeric`, even if
dichotomous. Methods we included and are conducted for feature computation
comprise autocorrelation analyses, (Bayesian) change-point analyses (`bcp`
[@erdman_bcp:_2007], `ecp` [@james_ecp_2015]), area determinations (`pracma`
[@pracma]), regression (multi-parametric non-linear & robust local & regression
models with segmented relationships: (`robustbase`
[@todorov_object-oriented_2009], `stats` [@stats], `segmented`
[@muggeo_interval_2017]) and hook effect detection (`PCRedux`
[@burdukiewicz_algorithms_2018]).

`encu()` is an extension of the `pcrfit_single()` function, where meta
information (e.g. detection chemistry, thermocycler) is included. It is suitable
for analysis of large data sets. All functions are protected against errors
(e.g. missing values). Both functions use, for instance, the following
descriptor-generating functions:

* `earlyreg()`, which calculates the features by regression analysis in the background region,
* `head2tailratio()`, which compares the ratio of head and tail,
* `hookregNL()` and `hookreg()`, which try to detect a hook effect [@burdukiewicz_algorithms_2018] in the amplification curve,
* `mblrr()` which performs a local robust regression analysis,
* `winklR()`, which is a function to calculate the angle based on the first and the second derivative of an amplification curve data from a quantitative PCR experiment and
* `autocorrelation_test()`, which is a function to test for autocorrelation of amplification curve data from a quantitative PCR experiment.

Auxiliary preprocessing and analysis functions of the package are:

* `armor()`, which is a helper function that catches errors and creates an output that can be used for further processing,
* `decision_mode()`, which calculates the frequency of classes in a dataset,
* `qPCR2fdata()`, which is helper function to convert amplification curve data to the *fdata* format for Hausdorff distance analysis [@fda.usc] and
* `performeR()`, which performs power analyses (e.g., sensitivity, specificity, Cohen's $\kappa$) for binary classification.

## Graphical User Interface:

`run_PCRedux()` is a graphical user interface (Figure \autoref{fig:fig_2})
based on `shiny` technology [@shiny] that outputs the features as a table. This
makes the results accessible to other programming languages and data analysis
tools.

![Graphical user interface for the analyses of qPCR data. A) The `run_PCRedux()`
GUI for analysis and tabular display can be used browsers or R environments that
support `ECMA Script` and `HTML`. In this example, the GUI was used in `RKWard`
(v.~0.7.2, Linux, Kubuntu 21.10, [@rodiger_rkward:_2012]). B) Depending on the
system used, information about the current analysis of errors can optionally be
output via the R console.\label{fig:fig_2}](fig_2.png)

## Datasets and Data Labeling

Unlabeled datasets are available in many laboratories. The function `tReem()`
was developed for the manual classification of amplification curves
(multicategory labeling). The first column must contain the qPCR cycles and all
subsequent columns must contain the amplification curves. Measures of similarity
between amplification curve shapes are the Pearson correlation coefficient or
the Hausdorff distance [@fda.usc].

The package contains labeled amplification curves (n = 14360; label: negative,
positive, ambiguous) from various qPCR instruments and detection methods, as
determined by the majority vote of an independent classification from four
experienced researchers (@PCRedux).

# Conclusion

Manual classification of amplification curves is time-consuming and error-prone,
especially when data sets are large. The majority of authors of other software
have so far been concerned with the determination of *Cq* value, amplification
efficiency and corresponding algorithms to quantify the amount of DNA. It is
difficult to analyze and classify amplification curves when they deviate
markedly from the sigmoidal shape or when their number becomes impractical for
manual analysis. The objectivity of the user can also be questioned because in a
manual classification (e.g., negative, positive) the result is usually
influenced by the subjective perception of the experimenter. An automatic system
for analyzing qPCR curves may objectify and generalize the decision process.

To address these problems, we developed `PCRedux`. With our package, 90
different descriptors of amplification curves can be calculated. The `PCRedux`
package allows an automatic computer-aided classification of amplification
curves based on descriptors derived from statistically and analytically
comprehensible algorithms. Especially for diagnostic applications, we see
advantages compared to neural networks. The analysis process is made faster
(omission of manual inspection), more objective and more reproducible. Examples
of application in the context of machine learning (classification positive vs.
negative) can be found in @PCRedux vignette. We hope that the examples will show
how the development of pipelines for machine learning and automatic
classification of amplification curves is possible.

# Acknowledgments

Grateful thanks belong to the `R` community.

# Funding

None

# References
