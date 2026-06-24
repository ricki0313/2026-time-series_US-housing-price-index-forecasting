# Time Series Analysis and Forecasting of the U.S. Housing Price Index
## Abstract
This project analyzes the U.S. Housing Price Index (HPI) using time series models to describe its historical patterns and forecast short-term trends. The analysis shows that HPI has a clear long-term upward trend, with faster growth after 2002, followed by a significant decline and higher volatility during the 2007–2010 financial crisis. After applying log transformation and first differencing, this project compares ARIMA, SARIMA, ARIMAX, and SARIMAX models using AICc, BIC, residual diagnostics, and out-of-sample forecast errors. The results show that HPI monthly growth has short-term persistence, a 12-period dependence structure, and that the crisis dummy has a significant negative effect on growth. Considering model fit, forecasting performance, and interpretability, SARIMAX(1,1,1)(1,0,0)[12] with crisis dummy is selected as the final model. The forecast results show that the model captures the gradual recovery of HPI after 2012, although it slightly underestimates the rapid increase in the later testing period.
---
## Dataset and Features
This project uses `monthly-housing.csv`, obtained from the Economics/Finance section of the Time Series Data Sets provided by Texas Tech University.

The dataset includes:
1. `date`: monthly observation date
2. `hpi`: U.S. Housing Price Index

The data period is from January 1991 to April 2013, with 268 valid monthly observations.  
This project mainly uses `hpi` as the target variable for analysis and forecasting, and the last 12 months are used as the testing set.
---
## Execution Order
Please run the files in the following order:
1. `01_data_cleaning.ipynb`
2. `02_eda.ipynb`
3. `03_stationary_inspection.ipynb`
4. `04_trans_diff.ipynb`
5. `05_candidate_model.R`
6. `06_model_diagnostics.R`
7. `07_forecasting.R`
---
For the complete analysis and results, please refer to `report.pdf` on the main page.