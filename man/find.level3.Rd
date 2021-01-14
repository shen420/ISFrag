\name{find.level3}
\alias{find.level3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finds Level 3 In-source Fragmentations
}
\description{
Creates a list of tables, where each table contains a parent feature and all its Level 1 in-source fragment features.
}
\usage{
find.level3(MS1directory, MS1.files, featureTable, type, peakCOR = 0.8, loss = 10, mz.tol = 0.02, rt.tol = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MS1directory}{
Directory containing .mzXML files in which MS1 scans were previously extracted from. This string is outputted by the generate.featuretable() function.
}
  \item{MS1.files}{
String vector containing names of .mzXML files in which MS1 scans were previously extracted from. This vector is outputted by the generate.featuretable() function.
}
  \item{featureTable}{
FeatureTable containg MS2 fragment information. This table is produced by the ms2.tofeaturetable() function.
}
  \item{type}{
For single sample analysis: type = "single", for multiple samples analysis: type = "multi".
}
  \item{peakCOR}{
MS1 EIC Peak-peak Pearson correlation threshold to determine if two features are co-eluting, and thus a Level 1 in-source fragment. Default value is 0.8.
}
  \item{loss}{
Mass of in-source fragment has to be "loss" dalton lower than that of its parent. Default value for loss is 10 daltons.
}
  \item{mz.tol}{
Mass range used to extract feature EIC. Default value is 0.02 daltons.
}
  \item{rt.tol}{
Retention time range used to extract feature EIC. Default value is 15 seconds.
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
function (MS1directory, MS1.files, featureTable, type, peakCOR = 0.8,
    loss = 10, mz.tol = 0.02, rt.tol = 15)
{
    peak_smooth <- function(x, level = 2) {
        n <- level
        if (length(x) < 2 * n) {
            return(x)
        }
        else if (length(unique(x)) == 1) {
            return(x)
        }
        else {
            y <- vector(length = length(x))
            for (i in 1:n) {
                y[i] <- sum(c((n - i + 2):(n + 1), n:1) * x[1:(i +
                  n)])/sum(c((n - i + 2):(n + 1), n:1))
            }
            for (i in (n + 1):(length(y) - n)) {
                y[i] <- sum(c(1:(n + 1), n:1) * x[(i - n):(i +
                  n)])/sum(c(1:(n + 1), n:1))
            }
            for (i in (length(y) - n + 1):length(y)) {
                y[i] <- sum(c(1:n, (n + 1):(n + i - length(x) +
                  1)) * x[(i - n):length(x)])/sum(c(1:n, (n +
                  1):(n + i - length(x) + 1)))
            }
            return(y)
        }
    }
    no_cores <- detectCores() - 2
    print("Using cores:")
    print(no_cores)
    registerDoParallel(no_cores)
    if (type == "single") {
        setwd(MS1directory)
        xraw <- xcmsRaw(MS1.files, profstep = 0)
        ISFtable <- foreach(i = (1:nrow(featureTable)), .packages = c("xcms",
            "MSnbase", "dplyr")) \%dopar\% {
            currFeature <- featureTable[i, ]
            minRT <- currFeature$rt - rt.tol
            maxRT <- currFeature$rt + rt.tol
            similarFeatures <- featureTable[featureTable$rt >
                minRT & featureTable$rt < maxRT, ]
            if (nrow(similarFeatures) == 0) {
                return(NA)
            }
            similarFeatures <- similarFeatures[similarFeatures$mz <=
                currFeature$mz[1] - loss, ]
            if (nrow(similarFeatures) == 0) {
                return(NA)
            }
            for (u in 1:nrow(similarFeatures)) {
                if ((sum(similarFeatures[u, 5:(5 + length(MS1.files) -
                  1)])/(length(MS1.files))) > (sum(currFeature[1,
                  5:(5 + length(MS1.files) - 1)])/(length(MS1.files)))) {
                  similarFeatures[u, ] <- NA
                }
            }
            similarFeatures <- similarFeatures[complete.cases(similarFeatures),
                ]
            if (nrow(similarFeatures) == 0) {
                return(NA)
            }
            curr.mass.lower.limit <- currFeature$mz[1] - mz.tol
            curr.mass.upper.limit <- currFeature$mz[1] + mz.tol
            curr.rt.lower.limit <- currFeature$rt[1] - rt.tol
            curr.rt.upper.limit <- currFeature$rt[1] + rt.tol
            if (curr.rt.lower.limit > tail(xraw@scantime, n = 1) |
                curr.rt.upper.limit > tail(xraw@scantime, n = 1)) {
                return(NA)
            }
            if (curr.rt.lower.limit < xraw@scantime[1] + 1) {
                curr.rt.lower.limit <- xraw@scantime[1] + 1
            }
            if (curr.rt.lower.limit < 1) {
                curr.rt.lower.limit <- 1
            }
            if (curr.rt.upper.limit > tail(xraw@scantime, n = 1)) {
                curr.rt.upper.limit <- tail(xraw@scantime, n = 1) -
                  1
            }
            if (curr.mass.lower.limit < xraw@mzrange[1]) {
                return(NA)
            }
            if (curr.mass.upper.limit > xraw@mzrange[2]) {
                return(NA)
            }
            mzRange <- as.double(cbind(curr.mass.lower.limit,
                curr.mass.upper.limit))
            RTRange <- as.integer(cbind(curr.rt.lower.limit,
                curr.rt.upper.limit))
            eeic <- rawEIC(xraw, mzrange = mzRange, rtrange = RTRange)
            currEIC <- eeic[["intensity"]]
            putativeISF <- data.frame(matrix(ncol = ncol(featureTable) +
                1, nrow = 0))
            ppcor <- 0
            putativeISF <- rbind(putativeISF, cbind(currFeature,
                ppcor))
            for (j in 1:nrow(similarFeatures)) {
                mass.lower.limit <- similarFeatures$mz[j] - mz.tol
                mass.upper.limit <- similarFeatures$mz[j] + mz.tol
                rt.lower.limit <- curr.rt.lower.limit
                rt.upper.limit <- curr.rt.upper.limit
                if (mass.lower.limit < xraw@mzrange[1])
                  (next)()
                if (mass.upper.limit > xraw@mzrange[2])
                  (next)()
                mzRange <- as.double(cbind(mass.lower.limit,
                  mass.upper.limit))
                RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
                eeic <- rawEIC(xraw, mzrange = mzRange, rtrange = RTRange)
                tmpEIC <- eeic[["intensity"]]
                currEIC <- currEIC[1:min(length(currEIC), length(tmpEIC))]
                tmpEIC <- tmpEIC[1:min(length(currEIC), length(tmpEIC))]
                if (is.na(cor(peak_smooth(currEIC), peak_smooth(tmpEIC))))
                  (next)()
                ppcor <- cor(peak_smooth(currEIC), peak_smooth(tmpEIC))
                if (ppcor >= peakCOR) {
                  putativeISF <- rbind(putativeISF, cbind(similarFeatures[j,
                    ], ppcor))
                }
            }
            if (nrow(putativeISF) == 1) {
                return(NA)
            }
            else {
                putativeISF$ISF_level[2:nrow(putativeISF)] <- "Level_3"
                putativeISF$ISF_level[1] == "Parent"
                return(putativeISF)
            }
        }
    }
    else {
        setwd(MS1directory)
        xraw <- list()
        for (w in 1:length(MS1.files)) {
            xraw[[w]] <- xcmsRaw(MS1.files[w], profstep = 0)
        }
        ISFtable <- foreach(i = (1:nrow(featureTable)), .packages = c("xcms",
            "MSnbase")) \%dopar\% {
            currFeature <- featureTable[i, ]
            minRT <- currFeature$rt - rt.tol
            maxRT <- currFeature$rt + rt.tol
            similarFeatures <- featureTable[featureTable$rt >
                minRT & featureTable$rt < maxRT, ]
            if (nrow(similarFeatures) == 0) {
                return(NA)
            }
            similarFeatures <- similarFeatures[similarFeatures$mz <=
                currFeature$mz[1] - loss, ]
            if (nrow(similarFeatures) == 0) {
                return(NA)
            }
            for (u in 1:nrow(similarFeatures)) {
                if ((sum(similarFeatures[u, 5:(5 + length(MS1.files) -
                  1)])/(length(MS1.files))) > (sum(currFeature[1,
                  5:(5 + length(MS1.files) - 1)])/(length(MS1.files)))) {
                  similarFeatures[u, ] <- NA
                }
            }
            similarFeatures <- similarFeatures[complete.cases(similarFeatures),
                ]
            if (nrow(similarFeatures) == 0) {
                return(NA)
            }
            putativeISF <- data.frame(matrix(ncol = ncol(featureTable) +
                1, nrow = 0))
            ppcor <- 0
            putativeISF <- rbind(putativeISF, cbind(currFeature,
                ppcor))
            for (j in 1:nrow(similarFeatures)) {
                ppcor <- 0
                count <- 0
                for (k in 1:length(MS1.files)) {
                  if (currFeature[1, 4 + k] == 0)
                    (next)()
                  if (similarFeatures[j, 4 + k] == 0)
                    (next)()
                  curr.mass.lower.limit <- currFeature$mz[1] -
                    mz.tol
                  curr.mass.upper.limit <- currFeature$mz[1] +
                    mz.tol
                  curr.rt.lower.limit <- currFeature$rt[1] -
                    rt.tol
                  curr.rt.upper.limit <- currFeature$rt[1] +
                    rt.tol
                  if (curr.rt.lower.limit > tail(xraw[[k]]@scantime,
                    n = 1) | curr.rt.upper.limit > tail(xraw[[k]]@scantime,
                    n = 1))
                    (next)()
                  if (curr.rt.lower.limit < xraw[[k]]@scantime[1] +
                    1) {
                    curr.rt.lower.limit <- xraw[[k]]@scantime[1] +
                      1
                  }
                  if (curr.rt.lower.limit < 1) {
                    curr.rt.lower.limit <- 1
                  }
                  if (curr.rt.upper.limit > tail(xraw[[k]]@scantime,
                    n = 1)) {
                    curr.rt.upper.limit <- tail(xraw[[k]]@scantime,
                      n = 1) - 1
                  }
                  if (curr.mass.lower.limit < xraw[[k]]@mzrange[1])
                    (next)()
                  if (curr.mass.upper.limit > xraw[[k]]@mzrange[2])
                    (next)()
                  mzRange <- as.double(cbind(curr.mass.lower.limit,
                    curr.mass.upper.limit))
                  RTRange <- as.integer(cbind(curr.rt.lower.limit,
                    curr.rt.upper.limit))
                  eeic <- rawEIC(xraw[[k]], mzrange = mzRange,
                    rtrange = RTRange)
                  currEIC <- eeic[["intensity"]]
                  mass.lower.limit <- similarFeatures$mz[j] -
                    mz.tol
                  mass.upper.limit <- similarFeatures$mz[j] +
                    mz.tol
                  rt.lower.limit <- curr.rt.lower.limit
                  rt.upper.limit <- curr.rt.upper.limit
                  if (mass.lower.limit < xraw[[k]]@mzrange[1])
                    (next)()
                  if (mass.upper.limit > xraw[[k]]@mzrange[2])
                    (next)()
                  mzRange <- as.double(cbind(mass.lower.limit,
                    mass.upper.limit))
                  RTRange <- as.integer(cbind(rt.lower.limit,
                    rt.upper.limit))
                  eeic <- rawEIC(xraw[[k]], mzrange = mzRange,
                    rtrange = RTRange)
                  tmpEIC <- eeic[["intensity"]]
                  currEIC <- currEIC[1:min(length(currEIC), length(tmpEIC))]
                  tmpEIC <- tmpEIC[1:min(length(currEIC), length(tmpEIC))]
                  if (is.na(cor(peak_smooth(currEIC), peak_smooth(tmpEIC))))
                    (next)()
                  ppcor <- ppcor + cor(peak_smooth(currEIC),
                    peak_smooth(tmpEIC))
                  count <- count + 1
                }
                if (count == 0)
                  (next)()
                ppcor <- ppcor/count
                if (ppcor >= peakCOR) {
                  putativeISF <- rbind(putativeISF, cbind(similarFeatures[j,
                    ], ppcor))
                }
            }
            if (nrow(putativeISF) == 1) {
                return(NA)
            }
            else {
                putativeISF$ISF_level[2:nrow(putativeISF)] <- "Level_3"
                putativeISF$ISF_level[1] <- "Parent"
                return(putativeISF)
            }
        }
    }
    stopImplicitCluster()
    ISFtable <- ISFtable[!is.na(ISFtable)]
    for (a in 1:length(ISFtable)) {
        names(ISFtable)[a] <- paste0(rownames(ISFtable[[a]][1,
            ]), "_", round(ISFtable[[a]]$mz[1], digits = 2),
            "_", round(ISFtable[[a]]$rt[1], digits = 0))
    }
    return(ISFtable)
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
