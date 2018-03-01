\name{PCRedux_datasets}
\alias{data_sample}
\alias{RAS002}
\alias{RAS002_decisions}

\title{The datasets implemented in PCRedux}

\description{
A compilation of datasets for method evaluation/comparison.
}

\usage{
data_sample
RAS002
RAS002_decisions
}

\details{
\bold{data_sample}\cr
Setup: Amplification curve data were analyzed with the encu() and the decision_modus() functions.\cr
Details:\cr
Data sets: batsch1, boggy, C126EG595, competimer, dil4reps94, guescini1, karlen1, lievens1, reps384, rutledge, testdat, vermeulen1, VIMCFX96_60, stepone_std.rdml, RAS002.rdml, RAS003.rdml, HCU32_aggR.csv, lc96_bACTXY.rdml.

\bold{RAS002}\cr
Setup: Amplification curve data of the RAS002.rdml data set.\cr
Details:\cr
Data sets: RAS002.rdml.

\bold{RAS002_decisions}\cr
Setup: Classes of the amplification curves from the RAS002.rdml data set.\cr
Details:\cr
Data sets: decision_res_RAS002.csv.
}

\author{
Stefan Roediger
}

\references{
Roediger, S., Burdukiewicz, M., Spiess, A.-N. & Blagodatskikh, K. Enabling reproducible real-time quantitative PCR research: the RDML package. \emph{Bioinformatics} (2017). doi:10.1093/bioinformatics/btx528\cr
Roediger, S., Burdukiewicz, M. & Schierack, P. chipPCR: an R package to pre-process raw data of amplification curves. \emph{Bioinformatics} 31, 2900--2902 (2015)\cr
Ritz, C. & Spiess, A.-N. qpcR: an R package for sigmoidal model selection in quantitative real-time polymerase chain reaction analysis. \emph{Bioinformatics} 24, 1549--1551 (2008).

}


\examples{
\dontrun{
## 'data_sample' dataset.
head(data_sample)

## 'RAS002.rdml' dataset as rda file.
head(RAS002)
head(RAS002_decisions)
}
}

\keyword{models}
