% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genes2.R
\name{plot_genes2}
\alias{plot_genes2}
\title{Plot gene locations for genomic interval}
\usage{
plot_genes2(genes,  scan1output , snpinfo, snp_lod, max_peak, haplo_peak_drop = NULL,
            minrow=4, padding=0.2,
            colors=c("black", "red3", "green4", "blue3", "darkorange2"),
            scale_pos=1, start_field="start", stop_field="stop",
            strand_field="strand", name_field="Name", ...)
}

\description{
Modified version of qtl2's \code{\link{=plot_genes}{plot_genes()}} function to allow plotting of the haplotype's QTL peak (assuming the region is set around the peak) and to easily identify which gene the SNP with the highest LOD score overlap. See \code{\link[=plot_snpasso2]{plot_snpasso2()}} description for more info.
}


\arguments{
\item{genes}{Data frame containing \code{start} and \code{stop} in
Mbp, \code{strand} (as \code{"-"}, \code{"+"}, or \code{NA}), and
\code{Name}.}

\item{scan1output}{Output of \code{\link[=scan1]{scan1()}} using
SNP probabilities derived by
\code{\link[=genoprob_to_snpprob]{genoprob_to_snpprob()}}.}

\item{snpinfo}{Data frame with SNP information with the following
columns (the last three are generally derived from with
\code{\link[=index_snps]{index_snps()}}):
\itemize{
\item \code{chr} - Character string or factor with chromosome
\item \code{pos} - Position (in same units as in the \code{"map"}
attribute in \code{genoprobs}.
\item \code{sdp} - Strain distribution pattern: an integer, between
1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
binary encoding indicates the founder genotypes
\item \code{snp} - Character string with SNP identifier (if
missing, the rownames are used).
\item \code{index} - Indices that indicate equivalent
groups of SNPs.
\item \code{intervals} - Indexes that indicate which marker
intervals the SNPs reside.
\item \code{on_map} - Indicate whether SNP coincides with a marker
in the \code{genoprobs}
}}

\item{snp_lod}{Same as 'scan1output'.The scan1output may have been modified somewhere so this extra parameter was created to obtain the original SNP's LOD scores.}

\item{max_peak}{Numeric value in Mbps of the location of the haplotype's QTL peak. This can be obtained from \code{\link[=max_scan1]{max_scan1()}} 'pos' column.}

\item{haplo_peak_drop}{A vector of length two containing the range of Mbps that have LOD scores within the QTL peak drop. This vector should be the output from \code{\link[=lod_drop_rng]{lod_drop_rng()}}.}

\item{minrow}{Minimum number of rows of genes in the plot}

\item{padding}{Proportion to pad with white space around the genes}

\item{colors}{Vectors of colors, used sequentially and then re-used.}

\item{scale_pos}{Factor by which to scale position (for example, to convert basepairs to Mbp)}

\item{start_field}{Character string with name of column containing the genes' start positions.}
\item{stop_field}{Character string with name of column containing the genes' stop positions.}

\item{strand_field}{Character string with name of column containing the genes' strands.}
\item{name_field}{Character string with name of column containing the genes' names.}

\item{...}{Optional arguments passed to \code{\link[graphics:plot]{graphics::plot()}}.}
}
\value{
None.
}
\description{
Plot gene locations for a genomic interval, as rectangles with gene
symbol (and arrow indicating strand/direction) below.
}
\section{Hidden graphics parameters}{

Graphics parameters can be passed via \code{...}. For
example, \code{xlim} to control the x-axis limits.
These are not included as formal
}
