# HighFreqDataEpi
Data and source for the manuscript "Leveraging high frequency data for epidemiological analysis"

# Leveraging High Frequency Data for Epidemiological Analysis
Contributors: Esther Li Wen Choo<sup>1</sup>, Guo Peihong<sup>1</sup>, Pranav Tewari<sup>1</sup>, Borame Sue Lee Dickens<sup>2</sup>, Jue Tao Lim<sup>1</sup>

# Affiliations
<sup>1</sup>Lee Kong Chian School of Medicine, Nanyang Technological University, Singapore  
<sup>2</sup>Saw Swee Hock School of Public Health, National University of Singapore, Singapore

# Motivation and Objectives
Understanding the key determinants of diseases is crucial for policymakers to formulate evidence-based public health policies. Epidemiological analysis facilitates our understanding of factors correlated with increased disease transmission. Current regression-based methods are limited by their inability to incorporate data of different frequencies. Hence, we explored the use of mixed data sampling (MIDAS) regression for modelling dengue cases in Singapore. MIDAS models can incorporate predictors which are of higher frequency than the response variable. We constructed Bayesian autoregressive, Bayesian LASSO, Bayesian adaptive LASSO and Bayesian group LASSO MIDAS models and fitted the models. Custom Markov chain Monte Carlo algorithms were developed to estimate the parameters of these models. We found positive associations between relative humidity, absolute humidity, air temperature, total precipitation and 1-week ahead dengue cases.  A 1g/m<sup>3</sup> rise in absolute humidity and 1% rise in relative humidity over 40 days was associated with a 4.33 (95% CrI 3.49, 4.83) and 2.22 (95% CrI 2.01, 2.58) increase in 1-week ahead cases respectively. Additionally, a 0.01 mm increase in total precipitation and 1&deg;C increase in air temperature was associated with a 2.21 (95% CrI 1.97, 2.45) and 6.63 (95% CrI 5.97, 7.78) increase in 1-week ahead cases respectively. Additionally, the MIDAS weights showed the varying contribution of predictor lags to dengue case counts. We demonstrated that incorporating high frequency predictors with MIDAS has practical value in public health by utilising high frequency information to understand factors associated with disease transmission. The use of MIDAS models can be expanded to other infectious diseases where high frequency predictors are available.
