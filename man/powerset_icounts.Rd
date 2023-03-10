% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{powerset_icounts}
\alias{powerset_icounts}
\title{Powerset intersection counts}
\usage{
powerset_icounts(df)
}
\arguments{
\item{df}{A 0/1 \link{data.frame} or \link{matrix}.
E.g. a (patient, omic) has a value of \code{1} if the specific patient has that
particular data modality, otherwise \code{0}.
Column names should be the name of the different omics in that case.}
}
\value{
a \link{tibble} with rows different omic combinations and columns:
\itemize{
\item \code{combo_name} => collapsed omic name (hyphen '-' is used for concatenation)
\item \code{omics} => list of omic names that make up the combo
\item \code{n_omics} => number of omics in the combo
\item \code{intersect_count} => number of patients who have all omics in the combo
}
}
\description{
Use on a 0/1 matrix, where e.g. rows are patients and columns
are different data modalities (omics).
This function will return for all possible combination of omics,
the number of patients that have information on \strong{ALL} these omics
(intersection counts).
}
\examples{
library(dplyr)

set.seed(42)
m = matrix(sample(x = c(0,1), size = 20, replace = TRUE), ncol = 5, nrow = 4)
colnames(m) = LETTERS[1:5]
pics = powerset_icounts(m)
pics \%>\%
  arrange(desc(intersect_count)) \%>\%
  filter(n_omics > 3)

}
