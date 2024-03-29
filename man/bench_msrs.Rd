% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures.R
\name{bench_msrs}
\alias{bench_msrs}
\title{Measures for Benchmarking Survival Models}
\usage{
bench_msrs()
}
\value{
a \link[data.table:data.table]{data.table} with the suggested survival
measures used for benchmarking purposes
}
\description{
A convenient function that returns a table with measures that
can be used for benchmarking survival learners.
}
\details{
The following measures are provided via \link{mlr3proba}:
\itemize{
\item Discrimination measures (input is ranking prediction - \code{crank})
\itemize{
\item \link[mlr3proba:mlr_measures_surv.cindex]{Harrell's C-index} (doesn't account
for censoring)
\item \link[mlr3proba:mlr_measures_surv.cindex]{Uno's C-index} (estimates censoring
distribution using Kaplan-Meier)
}
\item Calibration measures (input is survival predictions - \code{distr})
\itemize{
\item \link[mlr3proba:mlr_measures_surv.graf]{Integrated Brier Score} (with
\emph{proper = TRUE, method = 2}, estimates the censoring distribution using
Kaplan-Meier, is also a measure of discrimination)
\item \link[mlr3proba:mlr_measures_surv.rcll]{Right-Censored Log-Loss} (proper loss score)
\item \link[mlr3proba:mlr_measures_surv.dcalib]{D-Calibration} (with \emph{B = 10,
chisq = FALSE}, returning the measure and not the p-value)
}
}

The Integrated Brier Score and Right-Censored Log-Loss are also offered in
Explained Residual Variation (ERV) versions, where the resulting score can
be interpreted as the percentage increase in performance over the baseline
Kaplan-Meier (KM) model.
So the best possible value for ERV-standardized scores is 1 and negative
values correspond to worse performance than KM.

Any of the above measures \strong{can be used for model evaluation} (but depends
on several other factors if they \emph{should} be used, e.g. censoring
distribution, dataset size, proportion of events in the test set, etc.).
For model \strong{optimization/tuning}, we suggest the Right-Censored Log-Loss
for survival predictions and Uno's C-index for risk predictions
(e.g. linear predictors).
}
\examples{
library(mlr3proba)
ms = bench_msrs()

# Available measures ids
ms$id

# Get specific mlr3 measure objects
ms[id \%in\% c('ibrier','harrell_c')]$measure

}
