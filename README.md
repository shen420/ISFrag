ISFrag R Package User Manual
================
Sam Shen, Jian Guo, Tao Huan
24/03/2021

-   [Part 1: Introduction and
    Installation](#part-1-introduction-and-installation)
-   [Part 2: MS1 Feature Extraction](#part-2-ms1-feature-extraction)
    -   [2.1 XCMS Feature Extraction](#21-xcms-feature-extraction)
    -   [2.2 Additional Featuretable
        Input](#22-additional-featuretable-input)
-   [Part 3: MS2 Annotation](#part-3-ms2-annotation)
-   [Part 4: In-source Fragments
    Identification](#part-4-in-source-fragments-identification)
-   [Part 5: Results Export](#part-5-results-export)
    -   [5.1 Export ISF Relationship
        Tree](#51-export-isf-relationship-tree)
    -   [5.2 Export ISF Result
        Featuretable](#52-export-isf-result-featuretable)
-   [Part 6: Additional Details and
    Notes](#part-6-additional-details-and-notes)
-   [Part 2: QC](#part-2-qc)
-   [Part 3: Differential ATAC](#part-3-differential-atac)
-   [Part 4: GC bias](#part-4-gc-bias)
-   [Part 5: Differential analysis
    results](#part-5-differential-analysis-results)

# Part 1: Introduction and Installation

`ISFrag` is an R package for identifying and annotating in-source
fragments in LCMS metabolite featuretable. The package is written in the
language R and its source code is publicly available at
<https://github.com/HuanLab/ISFrag.git>

To install `ISFrag` package R version 4.0.0 or above is required, and we
recommend using RStudio to complete the installation and usage of
`ISFrag` by following the steps below:

``` r
# Install "devtools" package from CRAN if you do not already have it installed
install.packages("devtools")

# Load "devtools" package
library(devtools)

# Install "ISFrag" from Github using "devtools"
install_github("shen420/ISFrag")

# Load "ISFrag"
library(ISFrag)
```

# Part 2: MS1 Feature Extraction

## 2.1 XCMS Feature Extraction

## 2.2 Additional Featuretable Input

# Part 3: MS2 Annotation

# Part 4: In-source Fragments Identification

# Part 5: Results Export

## 5.1 Export ISF Relationship Tree

## 5.2 Export ISF Result Featuretable

# Part 6: Additional Details and Notes

*Now using `samples` make a plot showing the experimental design, with
time on the x axis, treatment on the y axis, and one plot on the left
and one on the right for the two replicates (e.g. using `facet_grid`).*

### `#?#` *Make the above plot. Each point should represent one of the samples. - 1 pt*

``` r
#here, if the point is there, it means such a sample exists, if absent it means that there is no such sample
```

*In this study, one of the things they were comparing was BRM014 to
DMSO. The drug BRM014 is dissolved in DMSO, so DMSO alone is the
appropriate control to gauge the effect of BRM014.*

### `#?#` *Can we compare BRM014 to DMSO across all time points? Why/why not? - 1 pt*

# Part 2: QC

*With most genomics data, it is important both that samples have
sufficient coverage, and that the samples have similar coverage. Either
case can lead to underpowered analysis, or misleading results. Calcualte
the read coverage for each sample. *

### `#?#` Make a plot with read coverage on the y-axis (total number of reads) and the samples on the x-axis. - 3 pt\*

``` r
# there are many ways you could do this; one of which is using the melt/cast functions from reshape
```

### `#?#` *Which sample has the most coverage? - 0.5 pt*

### `#?#` *Which sample has the least? - 0.5 pt*

### `#?#` *What is the % difference between the max and min (relative to the min)? - 0.5 pt*

*In cases where samples have vastly different coverage, you can
potentially down-sample the higher-coverage samples. Sometimes, throwing
out the data in this way can also introduce new problems, so we’re going
to stick with the data we have.*

*For this assignment, we will look only at BI\_protac vs control data. *

### `#?#` *Create a new data.frame containing only the BI\_protac and control samples - 1 pt*

### `#?#` *For this subset, calculate the counts per million reads (CPM) for each sample - 2 pt*

### `#?#` *Plot the kernel density estimate for CPM (x axis). 1 curve per sample, different colours per curve. - 1 pt*

### `#?#` *Plot the kernel density estimate for log(CPM+1) (x axis), coloured as before - 1 pt*

### `#?#` *Why do you think log-transforming is usually performed when looking at genomics data? What about adding 1 before log transforming? - 2 pt*

### `#?#` *Some regions have very large CPMs. Inspect the peaks for which CPM&gt;400. What do you notice about them? 3 pt*

*Normally, we would remove some of these regions before continuing (and
would redo the above steps). Since this is an assignment, we will
continue with the data as-is.*

*Often a good first step is to see if the data look good. One way to do
this is by seeing whether or not the signals in each sample correlate
with each other in ways you expect.*

### `#?#` *Calculate the pairwise correlations between log(CPM+1)s for the samples and plot them as a heatmap (samples x samples) - 3 pt*

### `#?#` *What do you expect the correlations between replicates to look like? Is that what you see? - 2 pt*

*It is common to exclude some regions from analysis. For instance, we
won’t be able to robustly identify those that are differential but have
low coverage even if they are truly differential, so there is no point
testing these. We will also remove mitochondrial regions, a common
contaminant of ATAC-seq data.*

### `#?#` *Filter your data, retaining only regions where the average counts per sample is greater than 10, and also remove mitochondrial regions - 3 pt*

### `#?#` *How many peaks did you have before? How many do you have now? - 1 pt*

# Part 3: Differential ATAC

*We want to know what regions are differentially accessible between
BI\_protac and the control.*

*Today, we’re going to use edgeR, which is designed for RNA-seq, but
works well on ATAC-seq as well. The user guide is here:*
<https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf>

### `#?#` *Make a count matrix called `countMatrix` for the BI\_protac and control samples, including only the peaks we retained above - 2 pt*

### `#?#` *Make an MA plot for allDEStatsPairedTreatControlvsProtac -2pt*

### `#?#` *Make an MA plot for allDEStatsPairedTime6vs24 - 1 pt*

*Now we’re going to test loess normalization instead.*

### `#?#` *Make the same two MA plots as before, but this time using the loess normalized analysis - 1 pt*

### `#?#` *What was the first normalization method? What changed in the MA plots? Which analysis do you think is more reliable and why? - 4 pt*

# Part 4: GC bias

*Next, we will look at potential GC bias in the data. We will again use
bioconductor *

### `#?#` *Extract the genomic DNA sequences for each peak using hg38 - 3 pt*

*Now we will see if there’s any relationship between peak CPM and GC
content for each of the samples.*

### `#?#` *Create scatter plots (one per sample, e.g. using facet\_wrap), including lines of best fit (GAM), where each plot shows GC content (x axis) vs CPM (y axis) for each peak (points) -2pt*

``` r
#please limit the y axis to between 0 and 50
```

### `#?#` *Repeat the above, but this time showing only the lines of best fit and all on the same plot - 2 pt*

### `#?#` *Given this result, predict whether we will see a significant relationship between GC content and logFC in our differential peak analysis (loess-normalized). Justify your prediction. Predicting “wrong” will not be penalized, as long as your justification is correct. Don’t retroactively change your answer. - 2 pt*

### `#?#` *Plot the relationship between GC and logFC for the loess-normalized ControlvsProtac analysis. Also include a line of best fit (blue) and y=0 (red) - 2 pt*

### `#?#` *Now plot the same thing for the NON loess-normalized ControlvsProtac analysis. - 1 pt*

### `#?#` *Was your prediction correct? Do you think we should also account for GC normalization in our differential ATAC analysis? Why/why not? - 3 pt*

*We will leave GC normalization as an optional exercise, and will not
actually do it here.*

# Part 5: Differential analysis results

### `#?#` *Suppose we perform the analyses above, redoing the differential analysis once more with GC normalization, and also considering that we tested loess and the default normalization methods. Did we P-hack? Why or why not? - 2 pt*

*Going forward, we will only use the initial analysis (**not loess
normalized**)*

### `#?#` *Now considering the two comparisons (6 vs 24 hours, and protac vs control). EdgeR performed a correction for MHT, but if we want to analyze the results from both comparisons, do we need to re-adjust to account for the fact that we tested two different hypothesis sets (time and treatment)? Why/not? - 2 pt*

### `#?#` *How many differential peaks did you find (FDR&lt;0.01). - 1 pt*

### `#?#` *Make a volcano plot of the allDEStatsPairedTreatControlvsProtac, with -log10(p-value) on the y axis and logFC on the x. Colour points that are significant at an FDR&lt;0.01. - 2 pt*

### `#?#` *Plot the logCPM (x axis) by -log10(Pvalue) (y axis), again colouring by FDR&lt;0.01. - 2 pt*

### `#?#` *Do you think our initial filtering on peaks with at least 10 reads on average per sample was a good choice? Why or why not?*

*At this point there are many other follow ups you can and would do for
a real differential analysis, but we leave these as optional exercises.
For example:* 1. Confirming that the differential peaks look correct
(e.g. CPM heatmap) 2. Confirming that peaks look differential on the
genome browser 3. Looking for motif enrichment 4. Performing a GREAT
analysis, including functional enrichment and assigning peaks to genes

*Knit your assignment as a github\_document and submit the resulting .md
and this .Rmd to your github, and complete the assignment submission on
Canvas. Make sure to include the graphs with your submission. *
