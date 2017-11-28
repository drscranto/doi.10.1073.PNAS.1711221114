# Predicting phenological shifts in a changing climate

## Abstract
Phenological shifts constitute one of the clearest manifestations of climate warming.  Advanced emergence is widely reported in high-latitude ectotherms, but a significant number of species exhibit delayed, or no change in, emergence.  Here we present a mechanistic theoretical framework that reconciles these disparate observations and predicts population-level phenological patterns based solely on data on temperature responses of the underlying life history traits.  Our model, parameterized with data from insects at different latitudes, shows that peak abundance occurs earlier in the year when warming increases the mean environmental temperature, but is delayed when warming increases the amplitude of seasonal fluctuations.  We find that warming does not necessarily lead to a longer activity period in high-latitude species because it elevates summer temperatures above the upper limit for reproduction and development.  Our findings both confirm and confound expectations for ectotherm species affected by climate warming: an increase in the mean temperature is more detrimental to low latitude species adapted to high mean temperatures and low-amplitude seasonal fluctuations; an increase in seasonal fluctuations is more detrimental to high-latitude species adapted to low mean temperatures and high-amplitude fluctuations. 

## Supplemental information: Data

Contains 3 data files from:

Lu et al 2010 https://www.jstage.jst.go.jp/article/aez/45/3/45_3_387/_article

Green Plant bug (Apolygus lucorum) - Temperate species
data_temp_R.csv

Amarasekare and Savage 2012 http://www.jstor.org/stable/10.1086/663677

Harlequin bug (Murganta histrionica) - Mediterranean species
data_med_R.csv

Dreyer and Baumgartner 1996 http://onlinelibrary.wiley.com/doi/10.1111/j.1570-7458.1996.tb00783.x/abstract

Pod sucking bug (Clavigralla shadabi) - Tropical species
data_trop_R.csv

## Supplemental information: Fitting thermal performance curves and seasonal temperature function (R code)

We estimate the parameters of the thermal responses of juvenile mortality, adult mortality, birth rate, and development rate for each of the three species. We use nls to fit the specific functional forms assumed in the text. For the temperate and Mediterranean species' responses of development to temperature, we lack data at the cold extreme to estimate two parameters: TL, the temperature at which enzyme activity is 50% of the maximum value, and the Arrhenius constant associated with development below this threshold. In order to make a best guess of these parameters we write a very simple sum of squares function and find parameters that minimize it locally. For TL, we use an initial value that is a linear approximation of the 50% threshold from the data.

## Supplemental information: Solving the population model DDEs (Python code)

The code in DDE_one_pop.py solves the syustem of equations for the population model for a given set of parameter values.


The function define_params defines all demographic parameters for three species, given as the only unout to the function spp which can be either 'trop', 'temp', or 'med'. It returns an array of all parameter values.

There are constant parameter values that never vary in our scenarios: a year has 365 days, climate change occurs over 100 years, the time step of the system of equations is one day, and sets one of the dde solver parameter values that affects behaviour around zero.

In each climate change scenario we must define the species affected (spp), the time over which to solve (max_years), the number of years to plot or save (keep_years), the change in mean temperature (delta_mean), the change in the amplitude of seasonal fluctuations (delta_ampl), the functional form of the response of competition strength to temperature (q_form), the process that is density dependent (dd), and solver parameters (tol, dde_dt, dde_hbsize).


The code relies on PyDDE and numpy modules. 


