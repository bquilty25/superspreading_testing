# Superspreading testing

This repository contains code and data for the analyses in [Disentangling the drivers of heterogeneity in SARS-CoV-2 transmission from data on viral load and daily contact rates](). The aim of this work was to combine data from contact surveys and viral load studies to estimate the secondary infection distribution of SARS-CoV-2 over the course of the pandemic in the UK in 2020 and the effectiveness of targeted testing strategies for reducing superspreading events.

## Installation

Clone/download this project onto your machine using the green button at the top right of this page.

The `pacman` R package, which is a package manager, is required to run the code. It can be installed with: 

```R
install.packages("pacman")
```

## Data
We used contact data from the [BBC Pandemic](https://github.com/adamkucharski/2020-cov-tracing/tree/62fe9be98e1d7ae7b49fd6fa0938f82970afb715/data) [1] and CoMix [2] contact surveys and [viral load trajectory parameter estimates](https://github.com/skissler/CtTrajectories_B117/tree/9a5b14eeb01d7c4b26eec80932d28eb3e9349ca1/output) from the literature [3] to simulate viral load trajectories. We converted viral load to infectiousness using laboratory data on the probability of culturing virus at different viral loads [4]. All data required for running the analyses is available in the [data](data) folder or the linked repositories.

## Running the code

The main analysis can be run by navigating to where the code was download and running:

```R
source("scripts/main.R")
```

The results can then be visualised by immediately following this with:

```
source("scripts/results.R")
```

## Output
The results plots are saved in the [results](results) folder.

## Built With

* [R version X.X.X (XXXX-XX-XX)](https://www.r-project.org/)

## Authors

* Billy Quilty: <billy.quilty@lshtm.ac.uk>

## References
1. Klepac, P. et al. Contacts in Context: Large-Scale Setting-Specific Social Mixing Matrices from the BBC Pandemic Project. medRxiv (2020). [http://doi.org/10.1101/2020.02.16.20023754](http://doi.org/10.1101/2020.02.16.20023754)

2. Gimma, A. et al. Changes in social contacts in England during the COVID-19 pandemic between March 2020 and March 2021 as measured by the CoMix survey: A repeated cross-sectional study. PLoS Med. 19, e1003907 (2022). [https://doi.org/10.1371/journal.pmed.1003907](https://doi.org/10.1371/journal.pmed.1003907)

3. Kissler, S. M. et al. Densely sampled viral trajectories suggest longer duration of acute infection with B.1.1.7 variant relative to non-B.1.1.7 SARS-CoV-2. medRxiv (2021). [https://doi.org/10.1101/2021.02.16.21251535](https://doi.org/10.1101/2021.02.16.21251535)

4. Pickering, S. et al. Comparative performance of SARS-CoV-2 lateral flow antigen tests and association with detection of infectious virus in clinical specimens: a single-centre laboratory evaluation study. Lancet Microbe 2, e461–e471 (2021). [https://doi.org/10.1016/S2666-5247(21)00143-9](https://doi.org/10.1016/S2666-5247(21)00143-9)
