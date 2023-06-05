# The Cox-Aalen Transformation Model
The class of Cox-Aalen transformation model with censored data

The transformation model (Zeng and Lin, 2006), provides a broad framework for modeling censored data. 
However, it has a limitation in that it assumes the baseline hazards function is the same across all individuals or groups. 
In reality, the baseline hazards can vary depending on the characteristics of the individuals or groups being studied.
To address this limitation, we propose a class of so-called Cox-Aalen transformation models that allow the baseline hazards function to depend on covariates. 
By incorporating this extension, The proposed model can capture the heterogeneity in the population more accurately. 
In addition, the proposed model includes both the transformation model and the Cox-Aalen model (Scheike and Zhang, 2002) as sepcial cases.
We provide R code which implements an Expectation-Solving (ES) algorithm for easy computations. 
The AMP dataset contains two harmonized randomized HIV Vaccine/Prevention trials (Corey et al., 2021).

References:

Corey, L., Gilbert, P. B., Juraska, M., Montefiori, D. C., Morris, L., Karuna, S. T., Edupuganti, S., Mgodi, N. M., deCamp, A. C., Rudnicki, E., et al. (2021). Two randomized trials of neutralizing antibodies to prevent HIV-1 acquisition. New England Journal of Medicine 384, 1003–1014.

Ning, X., Pan, Y., Sun, Y., and Gilbert, P. B. (2023). A semiparametric cox-aalen transformation model with censored data. Accepted at Biometrics.

Scheike, T. H. and Zhang, M.-j. (2002). An additive–multiplicative Cox–Aalen regression model. Scandinavian Journal of Statistics 29, 75–88.

Zeng, D. and Lin, D. (2006). Efficient estimation of semiparametric transformation models for counting processes. Biometrika 93, 627–640.

