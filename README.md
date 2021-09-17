# EVT-D-statistics-comparison

Python script to compare D statistics (Coles et al, 2001) among statistical models and plot the best model for each grid point. 
As an input data, it uses the negative log-likeloods (NLLH) of generalized pareto distribution models. NLLH are estimated using R ismev (https://CRAN.R-project.org/package=ismev) and extRemes packages (https://CRAN.R-project.org/package=extRemes) following Coles et al. (2001).

The script reproduces Figure 11 in Kim et al. (2021):

Kim, W. M., Blender, R., Sigl, M., Messmer, M., & Raible, C. C. (2021). Statistical characteristics of extreme daily precipitation during 1501 BCE–1849 CE in the Community Earth System Model. Climate of the Past Discussions, 1-38. https://doi.org/10.5194/cp-2021-61

Sample data used for the figure are provided as netCDF files, which contain NLLH of the stationary and non-stationary GPD models. 
Some post-processed datasets used in Kim et al. (2021) are available in Klima und Umweltphysik (KUP)-Universität Bern Nextcloud: 

https://cloud.climate.unibe.ch/s/HAXWy8BJMntonNi

Full datasets of CESM 1.2.2 1501BCE-1849CE simulations are available by request. Please contact to Woon Mi Kim (woonmi.kim@unibe.ch) or Prof. Dr. Christoph Raible (christoph.raible@unibe.ch)



To use the data or the script provided here and/or in KUP Nextcloud, please cite these papers:  

Kim, W. M., Blender, R., Sigl, M., Messmer, M., and Raible, C. C: Statistical characteristics of extreme daily precipitation during 1501 BCE–1849 CE in the Community Earth System Model. Climate of the Past Discussions, 1-38. https://doi.org/10.5194/cp-2021-61, 2021.

Gilleland,  E.  and  Katz,  R.  W.:  extRemes  2.0:  An  extreme  value  analysis  package  in  R,  Journal  of  Statistical  Software,  72,  1–39, https://doi.org/10.18637/jss.v072.i08, 2016.

Coles, S., Bawa, J., Trenner, L., and Dorazio, P.: An introduction to statistical modeling of extreme values, vol. 208, Springer, 2001.



