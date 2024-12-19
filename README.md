Bidirectional Impact of Periodontal Disease and Major Depression and Their Shared Genetic Mechanisms: Evidence from Clinical Observations and Genetic Analyses
==============================

This code is the details of the manuscript "Bidirectional Impact of Periodontal Disease and Major Depression and Their Shared Genetic Mechanisms: Evidence from Clinical Observations and Genetic Analyses".


## Table of Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)

 
## Overview

Background
Epidemiological evidence suggests a potential association between periodontal disease (PD) and major depression (MD). However, the bidirectional association and underlying genetic basis between these conditions remain insufficiently explored. This study aims to elucidate the bidirectional impact of PD and MD, along with their shared genetic mechanisms.  
Methods
Propensity score matching (PSM) was employed using UK Biobank data to examine the bidirectional association between PD and MD. Subsequently, utilizing genome-wide association study summary statistic assessed the bidirectional causal association and genetic correlations. Utilizing multimodal analysis results, a genetic scoring system, the Genetic Association Score (GAS), was developed to evaluate shared risk genes between the two diseases.
Results
Clinical studies have revealed a bidirectional association between PD and MD (MD risk: Odds ratio [OR] 1.51, 95% confidence interval [CI], 1.46-1.53, P<0.001; PD risk: OR, 1.51, 95% CI, 1.45-1.57, P<0.001). The bidirectional causal association was verified through Mendelian randomization (MD risk: OR, 1.19, 95% CI, 1.05 -1.34, P=0.006; PD risk: OR, 1.20, 95% CI, 1.12 -1.29, P<0.001). Genetically, significant genetic correlations was observed (β, 0.341, SE, 0.052). The GAS identified several potential shared risk genes, including SNCA and PSEN2. Genes may play pivotal roles in immune responses, cellular signaling, and metabolic adaptations.
Conclusions
This study demonstrates the bidirectional genetic correlation between PD and MD, identified potential shared risk genes. These findings provide insights into the underlying mechanisms of their comorbidity and offer valuable implications for clinical and public health applications.


## System Requirements
### Hardware requirements
Requires only a standard computer with enough RAM to support the in-memory operations.

### Software requirements
### OS Requirements
Supported for *windows* and *Linux*. The package has been tested on the following systems:
+ Windows: Win 11
+ Linux: Ubuntu 20.04

### R Dependencies
Mainly depends on the R 4.2.2  . The versions of packages are, specifically:

```
car ≥ 3.1-1
dplyr ≥ 1.1.0
stringr ≥ 1.5.0
cmprsk ≥ 2.2-11
dplyr ≥ 1.1.0
foreign ≥ 0.8-84
ggplot2 ≥ 3.4.1
ggsci ≥ 2.9
ggrepel ≥ 0.9.3
lava ≥ 1.7.2.1
Matching ≥ 4.10-8
mice ≥ 3.15.0
devtools ≥ 2.4.5
pec ≥ 2022.5.4
poLCA ≥ 1.6.0.1
plyr ≥ 1.8.8
prodlim ≥ 2019.11.13
reshape2 ≥ 1.4.4
rms ≥ 6.5-0
riskRegression ≥ 2022.11.28
survey ≥ 4.1.1
scales ≥ 1.2.1
survminer ≥ 0.4.9
survival ≥ 3.5-0
splines ≥ 4.1.2
timeROC ≥ 0.4


```

## Installation Guide

Installation via cloned repository:

```
$ git clone https://github.com/leescu/GAS.git
$ cd PD-HPL
```


