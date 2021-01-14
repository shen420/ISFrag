\name{generate.featuretable}
\alias{generate.featuretable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Table of Features
}
\description{
Generate feature table from raw .mzXML files using XCMS's CentWave peak detection algorithm. Perform feature alignment for multi-file analysis.
}
\usage{
generate.featuretable(MS1directory, MS2directory, type, CAMERA = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MS1directory}{
Directory of folder containing all .mzXML files containing MS1 scans. This folder must contain only relevant .mzXML files of either DDA, DIA, or Fullscan sample(s).
}
  \item{MS2directory}{
Directory of folder containing all .mzXML files containing MS2 scans. This folder must contain only relevant .mzXML files of DDA sample(s).
}
  \item{type}{
For single sample analysis: type = "single", for multiple samples analysis: type = "multi".
}
  \item{ppm}{
numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition. Default = 10.
}
  \item{peakwidth}{
numeric(2) with the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds. Default = c(10,120).
}
  \item{mzdiff}{
numeric(1) representing the minimum difference in m/z dimension required for peaks with overlapping retention times; can be negative to allow overlap. During peak post-processing, peaks defined to be overlapping are reduced to the one peak with the largest signal. Default = 0.01.
}
  \item{snthresh}{
numeric(1) defining the signal to noise ratio cutoff. Default = 6.
}
  \item{integrate}{
Integration method. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact. Default = 1.
}
  \item{prefilter}{
numeric(2): c(k, I) specifying the prefilter step for the first analysis step (ROI detection). Mass traces are only retained if they contain at least k peaks with intensity >= I. Default = c(3,100).
}
  \item{noise}{
numeric(1) allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection). Default = 100.
}
  \item{bw}{
Sample alignment: bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram. Default = 5.
}
  \item{mzwid}{
Sample alignment: width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples. Default = 0.015.
}
  \item{max}{
Sample alignment: maximum number of groups to identify in a single m/z slice. Default = 100.
}
  \item{CAMERA}{
Whether to perform CAMERA analysis to identify adducts and isotopes. "T" or "F".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (MS1directory, MS2directory, type, ppm = 10, peakwidth = c(10,
    120), mzdiff = 0.01, snthresh = 6, integrate = 1, prefilter = c(3,
    100), noise = 100, bw = 5, mzwid = 0.015, max = 100, CAMERA = F)
{
    cwp <- CentWaveParam(ppm = ppm, peakwidth = peakwidth, mzdiff = mzdiff,
        snthresh = snthresh, integrate = integrate, prefilter = prefilter,
        noise = noise)
    if (type == "single") {
        setwd(MS1directory)
        MS1.files <- list.files(pattern = ".mzXML")
        MS1data <- readMSData(MS1.files, msLevel. = 1, mode = "onDisk")
        MS1data <- findChromPeaks(MS1data, param = cwp)
        data_filtered <- filterMsLevel(MS1data, msLevel = 1L)
        xset <- as(data_filtered, "xcmsSet")
        setwd(MS2directory)
        MS2.files <- list.files(pattern = ".mzXML")
        MS2data <- readMSData(MS2.files, mode = "onDisk")
        featureTable <- as.data.frame(xset@peaks)
        featureTable <- featureTable[order(featureTable$mz, decreasing = TRUE),
            ]
        rownames(featureTable) <- paste("F", 1:nrow(featureTable),
            sep = "")
        if (CAMERA) {
            xsa <- xsAnnotate(xset)
            anF <- groupFWHM(xsa, perfwhm = 0.6)
            anI <- findIsotopes(anF, mzabs = 0.01)
            anIC <- groupCorr(anI, cor_eic_th = 0.75)
            anFA <- findAdducts(anIC, polarity = "positive")
            peaklist <- getPeaklist(anFA)
            peaklist <- peaklist[order(peaklist$mz), ]
            featureTable <- cbind(featureTable, peaklist$isotopes)
            colnames(featureTable)[ncol(featureTable)] <- "Isotopes"
            featureTable <- cbind(featureTable, peaklist$adduct)
            colnames(featureTable)[ncol(featureTable)] <- "Adduct"
            featureTable <- cbind(featureTable, as.numeric(peaklist$pcgroup))
            colnames(featureTable)[ncol(featureTable)] <- "pcgroup"
        }
    }
    else {
        setwd(MS1directory)
        MS1.files <- list.files(pattern = ".mzXML")
        MS1data <- readMSData(MS1.files, msLevel. = 1, mode = "onDisk")
        MS1data <- findChromPeaks(MS1data, param = cwp, SnowParam())
        data_filtered <- filterMsLevel(MS1data, msLevel = 1L)
        xset <- as(data_filtered, "xcmsSet")
        setwd(MS2directory)
        MS2.files <- list.files(pattern = ".mzXML")
        MS2data <- readMSData(MS2.files, mode = "onDisk")
        xset <- group(xset, bw = bw, minfrac = 0.5, mzwid = mzwid,
            minsamp = 1, max = max)
        xset <- retcor(xset, method = "obiwarp", profStep = 1)
        xset <- group(xset, bw = bw, minfrac = 0.5, mzwid = mzwid,
            minsamp = 1, max = max)
        xset <- fillPeaks(xset)
        XCMt <- data.frame(xset@groups)
        xcmI <- groupval(xset, value = "maxo")
        featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin,
            XCMt$rtmax, xcmI)
        colnames(featureTable)[1:4] <- c("mz", "rt", "rtmin",
            "rtmax")
        featureTable <- as.data.frame(featureTable)
        featureTable <- featureTable[order(featureTable$mz, decreasing = TRUE),
            ]
        rownames(featureTable) <- paste("F", 1:nrow(featureTable),
            sep = "")
        if (CAMERA) {
            xsa <- xsAnnotate(xset)
            anF <- groupFWHM(xsa, perfwhm = 0.6)
            anI <- findIsotopes(anF, mzabs = 0.01)
            anIC <- groupCorr(anI, cor_eic_th = 0.75)
            anFA <- findAdducts(anIC, polarity = "positive")
            peaklist <- getPeaklist(anFA)
            peaklist <- peaklist[order(peaklist$mz), ]
            featureTable <- cbind(featureTable, peaklist$isotopes)
            colnames(featureTable)[ncol(featureTable)] <- "Isotopes"
            featureTable <- cbind(featureTable, peaklist$adduct)
            colnames(featureTable)[ncol(featureTable)] <- "Adduct"
            featureTable <- cbind(featureTable, as.numeric(peaklist$pcgroup))
            colnames(featureTable)[ncol(featureTable)] <- "pcgroup"
        }
    }
    MS1.files <<- MS1.files
    MS2.files <<- MS2.files
    MS1data <<- MS1data
    MS2data <<- MS2data
    return(featureTable)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
