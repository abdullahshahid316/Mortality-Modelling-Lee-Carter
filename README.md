# Mortality-Modelling-Lee-Carter
This README document provides an overview of the ACST 3059 Assignment, focusing on mortality analysis and projection using statistical modeling in R.

Overview

The project involves using statistical methods and models to analyze mortality data for Australia, using the Human Mortality Database (HMD). The analysis includes extracting and processing the data, applying parametric and non-parametric curve fitting, and comparing models for the best fit.

Required Packages

demography
forecast
splines

Data Analysis

Exploratory data analysis with summary statistics and structure examination.
Mortality plot generation for the year 2019, including total, male, and female mortality.
Model Implementation

Parametric curve fitting using spline models.
Training, validation, and test data preparation.
Mean Squared Error (MSE) calculation for model comparison.
Implementation of smoothing spline models with hyperparameter tuning.
Model Comparison and Evaluation

Application of the Lee-Carter model for mortality rate projection.
Residual analysis and graphical representation of model fit.
Comparison of model projections up to 2050.
Results

The project includes visual comparisons of mortality rate projections for 2030, 2040, and 2050 using both the Lee-Carter model and smoothing spline models.
Mean Squared Error (MSE) calculations for model performance.
Graphs illustrating the projections and the fit of the models.
Conclusion

The assignment provides a comprehensive analysis and comparison of different statistical models for mortality rate projection. It demonstrates the application of statistical techniques in demography and actuarial science.

Instructions for Use

Install the required R packages.
Use the provided functions and scripts to extract and analyze the data.
Modify parameters and inputs as needed for specific analysis or data years.
