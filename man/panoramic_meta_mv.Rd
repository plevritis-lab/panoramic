\name{panoramic_meta_mv}
\alias{panoramic_meta_mv}
\title{Multilevel Random-Effects Meta-Analysis of PANORAMIC Features}
\usage{
panoramic_meta_mv(
  se,
  patient_col,
  group_col = NULL,
  sample_col = "sample",
  tau_structure = c("patient", "patient_sample"),
  method = "REML",
  group_tau2 = c("none", "separate"),
  warn_sigma2 = TRUE,
  vi_floor = NULL,
  sigma2_tol = 1e-06,
  sigma2_rel = 1e-04,
  control = NULL,
  sparse = FALSE,
  BPPARAM = BiocParallel::SerialParam()
)
}
\arguments{
\item{se}{A \code{SummarizedExperiment} with assays \code{yi} and \code{vi}.}

\item{patient_col}{Column in \code{colData(se)} identifying patients.}

\item{group_col}{Optional grouping column. When provided, pooled means are estimated with \code{mods = ~ 0 + group}. For two groups, contrast is \code{beta_diff = group2 - group1}; group order is factor levels (if factor) or lexicographic order of labels.}

\item{sample_col}{Optional sample-id column; if \code{NULL}, \code{colnames(se)} are used.}

\item{tau_structure}{Random-effects structure: \code{"patient"} or \code{"patient_sample"}.}

\item{method}{Estimator passed to \code{metafor::rma.mv()} (for example \code{"REML"}).}

\item{group_tau2}{If \code{"separate"}, additionally fit per-group multilevel models to estimate group-specific heterogeneity components.}

\item{warn_sigma2}{If \code{TRUE}, warn when many features have near-zero variance components.}

\item{vi_floor}{Optional handling for non-positive \code{vi}: \code{"median"}, \code{"group_median"}, or \code{"none"}.}

\item{sigma2_tol}{Absolute tolerance for near-zero \code{sigma2} diagnostics.}

\item{sigma2_rel}{Relative tolerance (fraction of median \code{vi}) for near-zero \code{sigma2} diagnostics.}

\item{control}{Optional named optimizer control list passed to \code{metafor::rma.mv()}.}

\item{sparse}{Logical forwarded to \code{metafor::rma.mv(sparse = ...)}.}

\item{BPPARAM}{BiocParallel backend used across features.}
}
\value{
The input \code{SummarizedExperiment} with meta-analysis columns appended to
\code{rowData(se)}. For grouped analyses, columns are group-prefixed
(\code{*_mu_hat}, \code{*_se_mu}, \code{*_p_mu}, \code{*_k}). For exactly two
groups, contrast columns \code{beta_diff}, \code{se_diff}, \code{z_diff},
\code{p_diff}, \code{fdr_diff} are added. Metadata are stored in
\code{metadata(se)$panoramic$meta_mv}.
}
\description{
Feature-wise multilevel random-effects pooling for PANORAMIC spatial statistics
using \code{metafor::rma.mv()}, accounting for clustering of samples within
patients.
}
\examples{
yi <- matrix(c(0.20, 0.15, 0.10, 0.05), nrow = 1)
vi <- matrix(c(0.04, 0.05, 0.06, 0.05), nrow = 1)
colnames(yi) <- colnames(vi) <- paste0("s", seq_len(4L))

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(yi = yi, vi = vi),
  rowData = S4Vectors::DataFrame(ct1 = "A", ct2 = "B", radius_um = 10),
  colData = S4Vectors::DataFrame(
    sample = paste0("s", seq_len(4L)),
    patient = paste0("p", seq_len(4L)),
    group = c("control", "control", "case", "case")
  )
)

se_mv <- panoramic_meta_mv(
  se = se,
  patient_col = "patient",
  group_col = "group",
  sample_col = "sample",
  BPPARAM = BiocParallel::SerialParam()
)

SummarizedExperiment::rowData(se_mv)$beta_diff
}
