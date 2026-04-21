\name{plot_forest}
\alias{plot_forest}
\title{Forest Plot for PANORAMIC Spatial Colocalization}
\usage{
plot_forest(
  se_meta,
  ct1,
  ct2,
  radius_um = NULL,
  group_col = "group",
  group_colors = NULL,
  text_size = 2,
  header_text_size = 2,
  title_text_size = 8,
  base_size = 5,
  show_est_se = TRUE,
  show_ci = TRUE
)
}
\arguments{
\item{se_meta}{SummarizedExperiment from \code{panoramic_meta_mv()} containing \code{yi}/\code{vi} assays and pooled columns in \code{rowData}.}

\item{ct1, ct2}{Cell-type labels selecting one feature.}

\item{radius_um}{Radius (microns). If \code{NULL}, the first available radius for the pair is used.}

\item{group_col}{Group column in \code{colData(se_meta)}.}

\item{group_colors}{Optional named color vector for groups.}

\item{text_size}{Text size for row annotations.}

\item{header_text_size}{Text size for table-like column headers.}

\item{title_text_size}{Title text size.}

\item{base_size}{Base size passed to \code{theme_classic()}.}

\item{show_est_se}{Show the \code{Est (SE)} text column.}

\item{show_ci}{Show the \code{95\% CI} text column.}
}
\value{
A \code{ggplot} forest-plot object.
}
\description{
Draw individual-sample estimates and pooled group estimates for one
cell-type pair at one radius.
}
\details{
This helper currently supports exactly two groups in \code{group_col}.
}
\examples{
yi <- matrix(c(0.12, 0.18, 0.35, 0.30), nrow = 1)
vi <- matrix(c(0.03, 0.04, 0.05, 0.05), nrow = 1)
colnames(yi) <- colnames(vi) <- paste0("s", seq_len(4L))

se_meta <- SummarizedExperiment::SummarizedExperiment(
  assays = list(yi = yi, vi = vi),
  rowData = S4Vectors::DataFrame(
    ct1 = "A",
    ct2 = "B",
    radius_um = 10,
    control_mu_hat = 0.15,
    control_se_mu = 0.10,
    case_mu_hat = 0.32,
    case_se_mu = 0.11
  ),
  colData = S4Vectors::DataFrame(
    group = c("control", "control", "case", "case")
  )
)

p <- plot_forest(
  se_meta,
  ct1 = "A",
  ct2 = "B",
  radius_um = 10,
  group_col = "group"
)

p
}
