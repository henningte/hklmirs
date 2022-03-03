
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hklmirs

This repository contains the data and code for our two manuscripts (in
preparation):

> Henning Teickner and Klaus-Holger Knorr (in preparation): Improving
> Models to Predict Holocellulose and Klason Lignin Contents for Peat
> Soil Organic Matter with Mid Infrared Spectra.

> Henning Teickner and Klaus-Holger Knorr (in preparation): Predicting
> Absolute Holocellulose and Klason Lignin Contents for Peat Remains
> Challenging.

### How to cite

Please cite this compendium as:

> Henning Teickner and Klaus-Holger Knorr, (2022). *Compendium of R code
> and data for “Improving Models to Predict Holocellulose and Klason
> Lignin Contents for Peat Soil Organic Matter with Mid Infrared
> Spectra” and “Predicting Absolute Holocellulose and Klason Lignin
> Contents for Peat Remains Challenging”*. Accessed 03 Mar 2022. Online
> at <https://github.com/henningte/hklmirs/>

## Contents

The **analysis** directory contains:

  - [:file\_folder: paper](/analysis/paper): R Markdown source documents
    needed to reproduce the manuscript, including figures and tables.
    The main script is
    [001-paper-main.Rmd](analysis/paper/001-paper-main.Rmd). This script
    produces both manuscripts and the corresponding supplementary
    information. Additional scripts are:
    
      - [002-paper-m-original-models.Rmd](analysis/paper/002-paper-m-original-models.Rmd):
        Computes the original models used in Hodgkins et al. (2018) and
        models with the same model structure, but as Bayesian models.  
      - [003-paper-m-gaussian-beta.Rmd](analysis/paper/003-paper-m-gaussian-beta.Rmd):
        Computes models assuming a Beta distribution for holocellulose
        and Klason lignin contents and compares them to the original
        models.
      - [004-paper-m-reduce-underfitting.Rmd](analysis/paper/004-paper-m-reduce-underfitting.Rmd):
        Extents the Beta regression models by including additional
        variables (additional peaks) or using a different approach
        (using measured spectral intensities of binned spectra instead
        of extracted peaks), and validates these models using LOO-CV.
      - [005-paper-m-minerals.Rmd](analysis/paper/005-paper-m-minerals.Rmd):
        Uses the models from `003-paper-m-gaussian-beta.Rmd` to test how
        accurate a model for holocellulose content is which is also
        calibrated on training samples with higher mineral contents.
      - [006-paper-m-prediction-domain.Rmd](analysis/paper/006-paper-m-prediction-domain.Rmd):
        Analyzes the prediction domain (Wadoux et al. 2021) of the
        original models and the modified models and identifie under
        which conditions models extrapolate for peat and vegetation
        smaples from Hodgkins et al. (2018).
      - [007-paper-m-prediction-differences.Rmd](analysis/paper/007-paper-m-prediction-differences.Rmd):
        Compares predictions for the training data and the peat and
        vegetation data from Hodgkins et al. (2018) for the original
        models from Hodgkins et al. (2018) and the modified models from
        `004-paper-m-reduce-underfitting.Rmd`.
      - [008-paper-supplementary.Rmd](analysis/paper/008-paper-supplementary.Rmd):
        Computes supplementary analyses and figures for the first
        manuscript.
      - [001-reply-main.Rmd](analysis/paper/001-reply-main.Rmd): This is
        the main script for manuscript 2. It is run from within
        `001-paper-main.Rmd` and produces the supplementary information
        for manuscript 2.
      - [002-reply-main.Rmd](analysis/paper/002-reply-main.Rmd): This
        script produces the document for manuscript 2. It is run from
        within `001-reply-main.Rmd`.

  - [:file\_folder: data](/analysis/data): Data used in the analysis.
    Note that raw data is not stored in [:file\_folder:
    raw\_data](/analysis/data/raw_data) (empty folder), but in
    [:file\_folder: /inst/extdata](/inst/extdata). [:file\_folder:
    derived\_data](/analysis/data/derived_data) contains derived data
    computed from the scripts. The raw data are derived from Hodgkins et
    al. (2018).

  - [:file\_folder: stan\_models](/analysis/stan_models): The Stan model
    used in `001-reply-main.Rmd`.

The other folders in this directory follow the standard naming scheme
and function of folders in R packages. There are the following
directories and files:

  - `README.md`/`README.Rmd`: Readme for the compendium.  
  - `DESCRIPTION`: The R package DESCRIPTION file for the compendium.  
  - `NAMESPACE`: The R package NAMESPACE file for the compendium.  
  - `LICENSE.md`: Details on the license for the code in the
    compendium.  
  - `CONTRIBUTING.md` and `CONDUCT.md`: Files with information on how to
    contribute to the compendium.  
  - `Dockerfile`: Dockerfile to build a Docker image for the
    compendium.  
  - `.Rbuildignore`, `.gitignore`, `.dockerignore`: Files to ignore
    during R package building, to ignore by Git, and to ignore while
    building a Docker image, respectively.  
  - `renv.lock`: renv lock file (Lists all R package dependencies and
    versions and can be used to restore the R package library using
    renv). `renv.lock` was created by running `renv::snapshot()` in the
    R package directory and it uses the information included in the
    `DESCRIPTION` file.  
  - `.Rprofile`: Code to run upon opening the R-project.  
  - `R`, `man`, `inst`, `data-raw`, `data`, `src`: Default folders for
    making the R package run.
  - Folder `inst/extdata`: Folder with the raw data used for the
    analyses. All files in this folder are derived from Hodgkins et al.
    (2018).

## How to run in your broswer or download and run locally

You can download the compendium as a zip from from this URL:
<https://github.com/henningte/hklmirs/>

Or you can install this compendium as an R package, hklmirs, from GitHub
with:

``` r
remotes::install_github("henningte/hklmirs")
```

## How to use

#### Reproduce the analyses

To reproduce the analyses for the paper, open the Rstudio project
included in this research compendium and run the Rmarkdown script in
`analysis/paper/001-paper-main.rmd`.

Running the whole script takes about 12 hours and occupies additional
disk space of \~2 Gb.

Alternatively, the [Dockerfile](Dockerfile) can be used to build a
Docker image from which all analyses can be reproduced. The Dockerfile
ensures that all required dependencies are installed (e.g. specific R
packages; this is managed using the R package
[renv](https://rstudio.github.io/renv/articles/renv.html)).

The [Dockerfile](Dockerfile) provides instructions how to build a Docker
image from the Dockerfile and how to run the image in a Docker
container. It occupies disk space of \~7 Gb.

When the Docker image runs in a container, go to `localhost:8787` in
your Browser. You will find an RStudio interface where you can log in
with username `rstudio` and password `hkl`. Here you can find the
Rmarkdown scripts (`hklmirs/analysis/paper/001-paper-main.rmd`) as
described above.

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse. See the sources section for licenses for
data derived from external sources and how to give credit to the
original author(s) and the source.

### Sources

All files in `inst/extdata` are derived from Hodgkins et al. (2018).
These data are licensed under the
[CC-BY 4.0](http://creativecommons.org/licenses/by/4.0/) license (see
<https://www.nature.com/articles/s41467-018-06050-2#rightslink>).

The format of this research compendium is inspired by Marwick,
Boettiger, and Mullen (2018) and was created with rrtools (Marwick
2019). The Rmarkdown template for the main article is from the rticles
package (Allaire et al. 2020).

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.

### References

<div id="refs" class="references">

<div id="ref-Allaire.2020">

Allaire, JJ, Yihui Xie, R Foundation, Hadley Wickham, Journal of
Statistical Software, Ramnath Vaidyanathan, Association for Computing
Machinery, et al. 2020. *Rticles: Article Formats for R Markdown*.
Manual.

</div>

<div id="ref-Hodgkins.2018">

Hodgkins, Suzanne B., Curtis J. Richardson, René Dommain, Hongjun Wang,
Paul H. Glaser, Brittany Verbeke, B. Rose Winkler, et al. 2018.
“Tropical Peatland Carbon Storage Linked to Global Latitudinal Trends
in Peat Recalcitrance.” *Nature Communications* 9 (1): 3640.
<https://doi.org/10.1038/s41467-018-06050-2>.

</div>

<div id="ref-Marwick.2019">

Marwick, Ben. 2019. “Rrtools: Creates a Reproducible Research
Compendium.”

</div>

<div id="ref-Marwick.2018">

Marwick, Ben, Carl Boettiger, and Lincoln Mullen. 2018. “Packaging Data
Analytical Work Reproducibly Using R (and Friends).” *The American
Statistician* 72 (1): 80–88.
<https://doi.org/10.1080/00031305.2017.1375986>.

</div>

<div id="ref-Wadoux.2021">

Wadoux, Alexandre M. J.-C., Brendan Malone, Budiman Minasny, Mario
Fajardo, and Alex B. McBratney. 2021. *Soil Spectral Inference with R:
Analysing Digital Soil Spectra Using the R Programming Environment*.
Progress in Soil Science. Cham: Springer International Publishing.
<https://doi.org/10.1007/978-3-030-64896-1>.

</div>

</div>
