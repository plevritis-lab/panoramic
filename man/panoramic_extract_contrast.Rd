\name{panoramic_extract_contrast}
\alias{panoramic_extract_contrast}
\title{Extract PANORAMIC Contrast, Meta, or Spatialstats Tables}
\usage{
panoramic_extract_contrast(
  se,
  feature_cols = c("ct1", "ct2", "radius_um"),
  what = c("contrast", "meta", "spatialstats"),
  drop_na = FALSE
)
}
\arguments{
\item{se}{A \code{SummarizedExperiment} from PANORAMIC workflow steps.}

\item{feature_cols}{Feature columns to include when available (for example \code{ct1}, \code{ct2}, \code{radius_um}).}

\item{what}{Which table to extract: \code{"contrast"}, \code{"meta"}, or \code{"spatialstats"}.}

\item{drop_na}{If \code{TRUE}, drop rows with missing extracted statistics.}
}
\value{
A \code{data.frame}. For \code{what = "contrast"}, the table includes
\code{beta_diff}, \code{se_diff}, \code{z_diff}, \code{p_diff}, and
\code{fdr_diff}, plus available feature columns.
}
\description{
Convenience helper for flattening PANORAMIC outputs into analysis-ready tables.
}
\examples{
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(dummy = matrix(0, nrow = 1, ncol = 1)),
  rowData = S4Vectors::DataFrame(
    ct1 = "A",
    ct2 = "B",
    radius_um = 10,
    beta_diff = 0.25,
    se_diff = 0.10,
    z_diff = 2.5,
    p_diff = 0.012,
    fdr_diff = 0.02
  )
)

panoramic_extract_contrast(se)
}
