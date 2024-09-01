
# Epyeregress
## Introduction
Epyeregress estimates the effective reproduction number from time series of reported case numbers of epidemics. It is a Python reimplementation of the method outlined by Jin et al. [<a href="https://doi.org/10.3390/v14071576">1</a>] to estimate the reproduction number R from infection data, available in the R package EpiRegress [<a href="https://github.com/ShihuiJin/EpiRegress">2</a>] .

 **Note:** The JAGS model is identical to that from  the original R implementation of [<a href="https://github.com/ShihuiJin/EpiRegress">Epiregress</a>], this is merely a python implementation using the same principles as the original code from Jin et al. [<a href="https://doi.org/10.3390/v14071576">1</a>].

## Dependencies
This project requires the following dependencies:
- **PyJAGS**: A Python interface to JAGS (Just Another Gibbs Sampler) used for Bayesian model sampling. 

  **Note:** `PyJAGS` is currently only supported on Linux. Ensure that you are working in a Linux environment to run this code.


## The main steps for estimation of the effective reproduction number are:

 1. **Data Preparation and Normalization**:
	 -  The input data: creates the covariate matrix (`accumulator_matrix`), from local cases, and imported cases, vaccination numbers, oxford scores and mobility numbers. 
	 - The serial interval distribution is computed based on user-specified mean and standard deviation, representing the time between successive infections in a transmission chain.
	 
 2. **Model Intialisation**:
	 -  Covariates are selected based on the user input or a default subset.
	-   The data is normalized to remove biases and to standardize the scale across different covariates.

 3. **Bayesian Model Specification**:
	 -  The model is defined using a combination of local and imported case incidences, with a focus on estimating Rt​ using a Bayesian framework.
	-   The likelihood of observing the number of local cases on any given day is modeled using a Negative Binomial distribution, which accounts for overdispersion in the data.

 4. **Sampling and Inference**:
	 -  The posterior distributions of model parameters, including Rt​, are sampled using Markov Chain Monte Carlo (MCMC) methods implemented in `PyJAGS`.
	 - The model's hyperparameters, such as the coefficients for the covariates, are estimated simultaneously, allowing for the inclusion of exogenous factors in explaining Rt.

 5. **Smoothing and Output**:
	 -  The raw Rt​ estimates are smoothed using a Gaussian filter to reduce variability and provide more interpretable results.
	-   The function returns summary statistics for Rt​, the covariate effects, and the mean case numbers.

 6. **Visualisation**:
	 -  The function includes utilities for plotting the covariate impacts on RtR_tRt​ and the time-varying estimates of Rt​, facilitating a clear understanding of the dynamics of transmission.
 

