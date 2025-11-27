Project Title: Advanced Landslide Susceptibility Modeling with FGAM and GAM

Developer: Vitor Recacho

Context: This computational framework, developed by Vitor Recacho, serves as the methodological core of the manuscript "On the use of rainfall time series for regional landslide prediction by means of functional regression" (Ahmed, Recacho, et al.).

Introduction: The framework presents a robust approach to Landslide Susceptibility Modeling (LSM) by rethinking how environmental triggers are processed. The core innovation of this code lies in its treatment of antecedent rainfall not as a scalar variable, but as a continuous function over time. By utilizing Functional Generalized Additive Models (FGAM), the script captures the complex, time-varying influence of precipitation on slope stability. These results are systematically compared against traditional Generalized Additive Models (GAM) to demonstrate the efficacy of functional regression in hazard prediction.

Methodological Highlights: The code executes the rigorous statistical workflow detailed in the manuscript, including:

Data Integration: Merging high-dimensional functional precipitation matrices with scalar topographical predictors.

Spatial Validation: A "Leave-One-Event-Out" cross-validation strategy to evaluate spatial transferability.

Robustness Testing: Randomized cross-validation loops to generate learning curves and test stability.

Visualization: Generation of the 3D mesh plots and performance curves presented in the study.
