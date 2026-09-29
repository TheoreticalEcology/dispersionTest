# Simulation design – dispersion tests

All data were simulated with `DHARMa::createData()` using one continuous predictor (`Environment1`). Unless stated otherwise: slope (`fixedEffects`) = 1 (the `createData` default), overdispersion = 0, and DHARMa residuals use `n = 250` simulations (the `simulateResiduals` default). Overdispersion levels "0–1" are `seq(0, 1, 0.1)`, which gives 11 levels. In GLMMs the sample size is the **total** number of observations (so observations per group = sample size / Ngroups). Random-intercept variance is 1 in every GLMM and 0 in every GLM.

## Main table

| Script | Aim | Model | Distribution | Ntrials (K) | Intercept | Slope | Sample size | Overdispersion | Ngroups | Tests evaluated | N combinations | N simulations per combination | Total datasets |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1_glm_pearsonChisq | Pearson statistic vs. χ² distribution (KS test) | GLM | Binomial | 10 | −3, −1.5, 0, 1.5, 3 | 1 | 10, 20, 50, 100, 200, 500, 1000, 10000 | 0 | – | KS test of Pearson χ² | 40 | 1,000 × 100 ⁽ᵃ⁾ | 4,000,000 |
| 1_glm_pearsonChisq | Pearson statistic vs. χ² distribution (KS test) | GLM | Poisson | – | −3, −1.5, 0, 1.5, 3 | 1 | 10, 20, 50, 100, 200, 500, 1000, 10000 | 0 | – | KS test of Pearson χ² | 40 | 1,000 × 100 ⁽ᵃ⁾ | 4,000,000 |
| 2_glm_type1 | Type I error | GLM | Binomial | 10 | −3, −1.5, 0, 1.5, 3 | 1 | 10, 20, 50, 100, 200, 500, 1000, 10000 | 0 | – | Pearson χ²; DHARMa sim-based; DHARMa refit (param. bootstrap) | 40 | **10,000** ⁽ᵇ⁾ | 400,000 |
| 2_glm_type1 | Type I error | GLM | Poisson | – | −3, −1.5, 0, 1.5, 3 | 1 | 10, 20, 50, 100, 200, 500, 1000, 10000 | 0 | – | Pearson χ²; DHARMa sim-based; DHARMa refit | 40 | **10,000** ⁽ᵇ⁾ | 400,000 |
| 3_glm_power | Power | GLM | Binomial | 10 | −3, −1.5, 0, 1.5, 3 | 1 | 10, 20, 50, 100, 200, 500, 1000, 10000 | 0–1 (11) | – | Pearson χ²; DHARMa sim-based; DHARMa refit | 440 | 1,000 | 440,000 |
| 3_glm_power | Power | GLM | Poisson | – | −3, −1.5, 0, 1.5, 3 | 1 | 10, 20, 50, 100, 200, 500, 1000, 10000 | 0–1 (11) | – | Pearson χ²; DHARMa sim-based; DHARMa refit | 440 | 1,000 | 440,000 |
| 3a_glm_dispersionStats | Effect of slope on dispersion statistics/power | GLM | Binomial | 10 | 0 | −2, −1, 0, 1, 2, 3 | 500 | 0–1 (11) | – | Pearson χ²; DHARMa sim-based; DHARMa refit | 66 | 1,000 | 66,000 |
| 3a_glm_dispersionStats | Effect of slope on dispersion statistics/power | GLM | Poisson | – | 0 | −2, −1, 0, 1, 2, 3 | 500 | 0–1 (11) | – | Pearson χ²; DHARMa sim-based; DHARMa refit | 66 | 1,000 | 66,000 |
| 3b_glmBin_ntrials | Effect of number of trials | GLM | Binomial | 5, 20 | 0 | 1 | 500 | 0–1 (11) | – | Pearson χ²; DHARMa sim-based; DHARMa refit | 22 | 1,000 | 22,000 |
| 4_glmm_pearsonChisq | Pearson χ² in GLMMs (two-sided and greater) | GLMM | Binomial | 10 | 0 | 1 | 1000 | 0–1 (11) | 10, 20, 50, 100 | Pearson χ² (two-sided, greater) | 44 | **10,000** ⁽ᵇ⁾ | 440,000 |
| 4_glmm_pearsonChisq | Pearson χ² in GLMMs (two-sided and greater) | GLMM | Poisson | – | 0 | 1 | 1000 | 0–1 (11) | 10, 20, 50, 100 | Pearson χ² (two-sided, greater) | 44 | **10,000** ⁽ᵇ⁾ | 440,000 |
| 5_glmm_power_bin_10 | Type I error and power | GLMM | Binomial | 10 | −3, −1.5, 0, 1.5, 3 | 1 | 50, 100, 200, 500, 1000 | 0–1 (11) | 10 | Pearson χ²; DHARMa uncond.; DHARMa cond.; DHARMa refit cond. | 275 | 1,000 | 275,000 |
| 5_glmm_power_bin_50 | Type I error and power | GLMM | Binomial | 10 | −3, −1.5, 0, 1.5, 3 | 1 | 100, 200, 500, 1000 | 0–1 (11) | 50 | same as above | 220 | 1,000 | 220,000 |
| 5_glmm_power_bin_100 | Type I error and power | GLMM | Binomial | 10 | −3, −1.5, 0, 1.5, 3 | 1 | 200, 500, 1000 | 0–1 (11) | 100 | same as above | 165 | 1,000 | 165,000 |
| 5_glmm_power_pois_10 | Type I error and power | GLMM | Poisson | – | −1.5, 0, 1.5, 3 | 1 | 50, 100, 200, 500, 1000 | 0–1 (11) | 10 | same as above | 220 | 1,000 | 220,000 |
| 5_glmm_power_pois_50 | Type I error and power | GLMM | Poisson | – | −3, −1.5, 0, 1.5, 3 | 1 | 100, 200, 500, 1000 | 0–1 (11) | 50 | same as above | 220 | 1,000 | 220,000 |
| 5_glmm_power_pois_100 | Type I error and power | GLMM | Poisson | – | −3, −1.5, 0, 1.5, 3 | 1 | 200, 500, 1000 | 0–1 (11) | 100 | same as above | 165 | 1,000 | 165,000 |
| 5_glmm_runtime_tests | Runtime | GLMM | Poisson | – | 0 | 1 | 1000 | 0.5 | 100 | Pearson χ²; DHARMa sim-based; DHARMa refit | 1 | 1,000 | 1,000 |
| 6_alternative_DHARMa | Alternative DHARMa test vs. number of DHARMa simulations (nSim = 5, 10, 50, 250, 1000) | GLM | Binomial | 10 | −1.5, 0, 1.5 | 1 | 100, 1000 | 0–1 (11) | – | Pearson χ²; DHARMa sim-based; alternative (approx. Pearson) | 330 | 1,000 | 330,000 |
| 6_alternative_DHARMa | Alternative DHARMa test vs. number of DHARMa simulations (nSim = 5, 10, 50, 250, 1000) | GLM | Poisson | – | −1.5, 0, 1.5 | 1 | 100, 1000 | 0–1 (11) | – | Pearson χ²; DHARMa sim-based; alternative (approx. Pearson) | 330 | 1,000 | 330,000 |
| 7_alternative_DHARma_glmm_power | Power of alternative DHARMa test (conditional) | GLMM | Binomial | 10 | −1.5, 0, 1.5 | 1 | 20, 50, 100, 200, 500, 1000 | 0–1 (11) | 10 | Alternative (approx. Pearson) | 198 | 1,000 | 198,000 |
| 7_alternative_DHARma_glmm_power | Power of alternative DHARMa test (conditional) | GLMM | Binomial | 10 | −1.5, 0, 1.5 | 1 | 100, 200, 500, 1000 | 0–1 (11) | 50 | Alternative (approx. Pearson) | 132 | 1,000 | 132,000 |
| 7_alternative_DHARma_glmm_power | Power of alternative DHARMa test (conditional) | GLMM | Binomial | 10 | −1.5, 0, 1.5 | 1 | 200, 500, 1000 | 0–1 (11) | 100 | Alternative (approx. Pearson) | 99 | 1,000 | 99,000 |
| 7_alternative_DHARma_glmm_power | Power of alternative DHARMa test (conditional) | GLMM | Poisson | – | −1.5, 0, 1.5 | 1 | 50, 100, 200, 500, 1000 | 0–1 (11) | 10 | Alternative (approx. Pearson) | 165 | 1,000 | 165,000 |
| 7_alternative_DHARma_glmm_power | Power of alternative DHARMa test (conditional) | GLMM | Poisson | – | −1.5, 0, 1.5 | 1 | 100, 200, 500, 1000 | 0–1 (11) | 50 | Alternative (approx. Pearson) | 132 | 1,000 | 132,000 |
| 7_alternative_DHARma_glmm_power | Power of alternative DHARMa test (conditional) | GLMM | Poisson | – | −1.5, 0, 1.5 | 1 | 200, 500, 1000 | 0–1 (11) | 100 | Alternative (approx. Pearson) | 99 | 1,000 | 99,000 |
| 8_glmm_approxDFpearson | Approximate df for Pearson χ² in GLMMs | GLMM | Poisson | – | −1.5, 0, 1.5 | 1 | 200, 500, 1000 | 0–1 (11) | 10, 50, 100 | Pearson χ² with df: naive, Satterthwaite, KR, KR2 | 297 | 1,000 | 297,000 |

⁽ᵃ⁾ In script 1, each parameter combination was simulated 1,000 times, and a single Kolmogorov–Smirnov test compared the 1,000 Pearson statistics with the χ² distribution with the corresponding residual degrees of freedom. This procedure was repeated 100 times (seeds 1–100), so the reported result is the proportion of significant KS tests out of 100 replicates (100,000 datasets per combination).

⁽ᵇ⁾ Scripts 2 and 4 used 10,000 simulations per combination. All other type I error and power simulations used 1,000.
