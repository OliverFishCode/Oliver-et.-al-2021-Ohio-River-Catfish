## File descriptions

This folder contains models used in analysis of catfish population dynamics and demographics. Models are run from within species-specific R scripts in the parent directory.

`population_analysis.stan` is a Stan model used to estimate parameters of von Bertalanffy growth functions, length-weight regressions, and gear-specific catch-curves. Other quantities (e.g. natural mortality) are derived from model posteriors in R scripts to reduce memory requirements. The ensemble of models runs simultaneously to provide naturally paired posterior draws that can be used to make posterior predictions of YPR.

`yprC.cpp` is a C++ file that contains calculations for estimating yield per recruit from Bayesian posteriors using the output from Stan model above, and following the approach of Dippold et al. (2016) for simulating slot regulations.
