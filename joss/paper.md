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
  - name: BTU Cottbus–Senftenberg, Faculty Environment and Natural Sciences, Senftenberg, Germany
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

qPCR (quantitative polymerase chain reaction) is indispensable in research,
diagnostics and forensics, because it provides quantitative information
about the amount of DNA in a sample [@pabinger_survey_2014]. The interpretation
of amplification curves is often difficult if the curve does not
does not follow a typical amplification curve {fig:fig_1}.

`PCRedux` is an `R` package [@stats] for feature extraction and classification
in the sense of explainable machine learning. It uses functions to compute
ninety Boolean and numerical descriptors from amplification curves using
classical statistical and analytical procedures.

`PCRedux` can also be used to determine *Cq* values and amplification
efficiencies (*E*) in high-throughput analysis.

In particular, open data are the cornerstones of science. Given the lack of data
sets with classified amplification curves, `PCRedux` includes functions for
aggregation, management and aggregate, manage and disseminate classified qPCR
datasets. (trichotomously classified by multiple qPCR users into _negative_,
_positive_ and and _unique_) from different qPCR instruments and detection
methods.

# Statement of need

Numerous data with sigmoidal ('S'-shaped) curves exist. An important use case
being amplification curve data from quantitative polymerase chain reactions
(qPCR). Among other use cases, qPCR is a widely used laboratory method for
precise detection and quantification of pathogens and gene expression analysis.
The latter contributed significantly to the understanding of physiological and
pathological processes in pharmacology, medicine and forensics.
[@pabinger_survey_2014; @kok_small_2018]

Available software packages provide workflows and criteria for processing qPCR
data (incl. preprocessing of raw data, fitting of non-linear models to raw
amplification curve data, calculation of threshold or second-derivative based
'cycle of quantification' (*Cq*) values, amplification efficiency (*E*), relative
gene expression analysis, normalization procedures and data management)
[@pabinger_survey_2014; @ruijter_evaluation_2013; @ruijter_efficiency_2021;
@ramakers_assumption-free_2003].

qPCR curves must fulfill quality criteria for analysis. Curves can be
categorized by the user on the basis of somewhat subjective criteria (e.g.,
sigmoidal shape, steepness, noisiness, "hook effect" presence)
[@burdukiewicz_algorithms_2018; @spiess_system-specific_2016;
@spiess_impact_2015; @hanschmann_looptag_2021_2]. Positive qPCR reactions
usually exhibit a sigmoidal shape, while negative amplification curves will have
a more flat and linear appearance curvatures (Figure \autoref{fig:fig_1}).

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
derivative, equals *Cq*), amptester_polygon (area under the curve) and cpDdiff (absolute
cycle distance between the maximum of the second and first derivative) were
calculated and plotted for the classes positive negative and ambiguous. Data
from htPCR dataset [@ritz_qpcr:_2008].\label{fig:fig_1}](fig_1.png)

The question arises as to how amplification curves can be objectively and
reproducibly assessed and automatically interpreted (e.g. as
positive/negative/unambiguous or low/high quality). In high-throughput
experiments, manual analysis is not readily feasible as errors are likely due to
exhaustion of the concentration capacity. Manual analyses are often not
reproducible because scientists in different laboratories use arbitrary
thresholds or make a subjective assessment of quality. Internal laboratory
guidelines seem to facilitate the amplification curves, but these cannot be
transferred to other laboratories. [@bustin_why_2010; @taylor_ultimate_2019;
@kim_experimenting_2018].

Automatically extracted features from amplification curves (e.g., *Cq*- and *E*
values, slopes, change points, features of local curve segments) could provide a
solid basis for machine learning classification. To date, there is no
open-source software that applies classical biostatistical methods for
explainable machine learning to amplification curves. Furthermore, there are no
class-labelled datasets that can be used in this context.

PCRedux" is the first open-source software that extracts 90 mathematical
descriptors (features) from raw amplification curves. Features are numerically
or analytically derived, quantifiable, informative properties of an
amplification curve in a scalar unit.

# Software engineering

`PCRedux` (v.\~1.1-2, [MIT license](https://mit-license.org/)) is an `R` package
(S3 class system). `R` was chosen because it provides comprehensive tools for reproducible
statistical and bioinformatics analyses [@gentleman_bioconductor:_2004;
@gentleman_statistical_2007; @rodiger_r_2015; @liu_r_2014;
@leeper_archiving_2014]. Unit tests using the `testthat` package [@testthat]
were used for software quality control of `PCRedux`.

## Functions

We divide amplification curves into regions of interest (ROI) for feature
calculation (vignette @PCRedux Figure 5). For example, the baseline phase is
located in the left tail region, while the plateau phase is located in the right
tail region. The exponential and linear phase is located between these two ROIs
(vignette @PCRedux Figure 5).

`PCRedux`'s algorithms, published by others and ourselves (`qpcR`
[@ritz_qpcr:_2008], `MBmca` [@rodiger_surface_2013], `chipPCR`
[@rodiger_chippcr:_2015]), were are adapted for qPCR analysis. The `PCRedux`
dependencies include packages for preprocessing of raw data (`chipPCR`,
`MBmca`), fitting of non-linear models to raw data, calculation of
quantification cycles and amplification efficiencies (`qpcR`). `pcrfit_single()`
is the main function for descriptor calculation of single amplification curves.
It generates a *data.frame* with 90 descriptors. All output values are of type
`numeric`, even if dichotomous. Methods we included for feature computation
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
(e.g. missing values). Both functions use, the following
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

Examples of application in the context of machine learning (classification
positive vs. negative) can be found in @PCRedux vignette.

## Graphical User Interface:

`run_PCRedux()` is a graphical user interface (figure \autoref{fig:fig_2}) based
on `shiny` technology [@shiny] that outputs the features as a table, making them
accessible to other programming languages and data analysis tools.

![Graphical user interface for the analyses of qPCR data. A) The `run_PCRedux()`
GUI for analysis and tabular display can be used browsers or R environments that
support `ECMA Script` and `HTML`. In this example, the GUI was used in `RKWard`
(v.~0.7.2, Linux, Kubuntu 21.10, [@rodiger_rkward:_2012]). B) Depending on the
system used, information about the current analysis of errors can optionally be
output via the R console.\label{fig:fig_2}](fig_2.png)

## Datasets and Data Labeling

The function `tReem()` facilitates the manual blinded random classification of
unlabeled amplification curves (multicategory labeling). Measures of similarity
between amplification curve shapes are the Pearson correlation coefficient or
the Hausdorff distance [@fda.usc].

`PCRedux` contains labeled amplification curves (n = 14360; label: negative,
positive, ambiguous) from various qPCR instruments and detection methods, as
determined by the majority vote of an independent classification from four
experienced researchers (@PCRedux).

# Conclusion

Manual classification of amplification curves is time-consuming and error-prone,
especially with large data sets. It is difficult to work amplification curves if
they deviate from the sigmoidal shape or their number is impractical for manual
analysis. With a manual classification, the result is influenced by subjective
perception. An automated system for analysing qPCR curves can objectify and
generalise the decision-making process.

Training neural networks is an option. However, the question here is what the
model considers relevant. Especially in diagnostic This is particularly
problematic for diagnostic applications. To avoid these 'black box' situations,
there are methods where each descriptor can be explained mathematically.

The `PCRedux` package enables a faster (elimination of manual inspection)
automatic computer-aided classification of amplification curves based on 90
descriptors using statistically and analytically traceable algorithms. PCRedux'
serves to improve the quality and reproducibility of the quality and
reproducibility of qPCR data analysis by systematic validation of input data.

# Acknowledgments

Grateful thanks belong to the `R` community.

# Funding

None

# References
