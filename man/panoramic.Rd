\name{panoramic}
\alias{panoramic}
\title{Run PANORAMIC End-To-End}
\usage{
panoramic(
  spe_list,
  design,
  cell_type = "cell_type",
  pairs = "auto",
  radii_um,
  stat = "local_comp_enrichment",
  nsim = 100,
  correction = NULL,
  min_cells = 5L,
  concavity = 50,
  window = c("concave", "convex", "rect"),
  group_col = "group",
  group_tau2 = c("none", "separate"),
  patient_col = NULL,
  sample_col = NULL,
  tau_structure = c("patient", "patient_sample"),
  method_mv = "REML",
  vi_floor = "group_median",
  seed = 123,
  boot = c("approx", "block"),
  tile_size = NULL,
  nx = NULL,
  ny = NULL,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
}
\arguments{
\item{spe_list}{Named list of SpatialExperiment objects (one per sample).}

\item{design}{data.frame with at least columns \code{sample} and \code{group}.}

\item{cell_type}{Character colData column containing cell-type labels.}

\item{pairs}{Either \code{"auto"} or a data.frame with columns \code{ct1}, \code{ct2}.}

\item{radii_um}{Numeric vector of radii (microns).}

\item{stat}{Spatial statistic passed to \code{panoramic_spatialstats()}.}

\item{nsim}{Integer bootstrap replicates per sample/pair.}

\item{correction}{Optional edge-correction method for spatstat-based statistics. Ignored when \code{stat = "local_comp_enrichment"}.}

\item{min_cells}{Minimum cells per type per sample during preparation.}

\item{concavity}{Concavity for concave-hull sample windows (larger values approach convex hulls).}

\item{window}{One of \code{"concave"}, \code{"convex"}, or \code{"rect"}.}

\item{group_col}{Group column in \code{colData} used for multilevel pooling and contrasts.}

\item{group_tau2}{If \code{"separate"}, fit additional per-group multilevel models to report group-specific heterogeneity components.}

\item{patient_col}{Optional patient-id column for multilevel pooling. If \code{NULL}, PANORAMIC uses \code{"patient"} from \code{design} when available, otherwise falls back to \code{sample_col}.}

\item{sample_col}{Optional sample-id column used in \code{panoramic_meta_mv()}. If \code{NULL}, PANORAMIC uses \code{"sample"} (aligned to \code{design$sample}).}

\item{tau_structure}{Random-effects structure in \code{panoramic_meta_mv()}: \code{"patient"} or \code{"patient_sample"}.}

\item{method_mv}{Estimator passed to \code{metafor::rma.mv()} (for example \code{"REML"}).}

\item{vi_floor}{Variance-flooring mode passed to \code{panoramic_meta_mv()}. Default is \code{"group_median"} (\code{"median"}, \code{"group_median"}, or \code{"none"}).}

\item{seed}{Optional random seed for reproducibility.}

\item{boot}{Bootstrap mode in \code{panoramic_spatialstats()}: \code{"approx"} or \code{"block"}.}

\item{tile_size}{Optional tile size for block bootstrap (same units as coordinates).}

\item{nx, ny}{Optional tile counts for block bootstrap (used when \code{boot = "block"}).}

\item{BPPARAM}{BiocParallel backend.}

\item{verbose}{Logical verbosity passed to bootstrap routines.}
}
\value{
A list with:
\itemize{
\item \code{prep}: output of \code{panoramic_prepare()}.
\item \code{stats}: output of \code{panoramic_spatialstats()}.
\item \code{pooled}: output of \code{panoramic_meta_mv()}.
\item \code{tables}: pre-flattened data.frames (spatialstats/meta/contrast).
}
}
\description{
\code{panoramic_analyze()} runs preparation, spatial statistics, and multilevel
pooling in one call. \code{panoramic()} is a convenience alias with identical
arguments and return structure.
}
\details{
\code{panoramic()} delegates directly to \code{panoramic_analyze()}.
}
\examples{
toy <- panoramic_simulate_dataset(seed = 1)
out <- panoramic(
  spe_list = toy$spe_list,
  design = toy$design,
  cell_type = "cell_type",
  radii_um = c(10, 20),
  nsim = 5,
  min_cells = 2,
  window = "rect",
  BPPARAM = BiocParallel::SerialParam()
)
names(out)
}
