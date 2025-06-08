# Extreme_Weather_Correlation

Hurricanes profoundly impact human populations and various industries such as the insurance sector. Accurate modelling of such extreme weather events is therefore critical for disaster preparedness and mitigation strategies, which can ultimately contribute to saving lives. This thesis focuses on the state of Florida, USA, with a particular emphasis on the city of Miami. We begin by exploring Conditional Autoregressive (CAR) models to analyse extreme weather data, specifically wind speed and total precipitation. While these variables are inherently correlated in real-world scenarios, standard CAR models do not accommodate such interdependencies. To address this limitation, we modify the CAR framework to incorporate the correlation structure between variables. The principal goal of this thesis is to extend this model further to include a third variable: surface pressure. The resulting model captures the complex interactions among extreme weather variables and provides informative insights at locations of interest such as Miami -- something the two-variable model failed to do well.

1. data_extraction.py: Contains the code to download the data we work with from the Climate Data Store.
2. EDA.ipynb: Spatially and temporally reduces the data followed by an analysis of the empirical correlation between variables during extreme weather events.
3. pre_icar.ipynb: Contains the necessary data preprocessing before fitting our model in Stan.
4. no_correlation_icar.R: Stan script for model assuming no correlation between spatial effects.
5. icar.R: Script for two variable CAR model that accounts for correlation.
6. three_variable.R: Script for implementing the model when adding an extra variable.
7. post_icar.ipynb: Taking our fitted parameters, then plotting and analysing the results. 
