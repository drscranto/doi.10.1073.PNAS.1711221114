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

## Supplemental information: Fitting thermal performance curves (R code)

We estimate the parameters of the thermal responses of juvenile mortality, adult mortality, birth rate, and development rate for each of the three species. We use nls to fit the specific functional forms assumed in the text. For the temperate and Mediterranean species' responses of development to temperature, we lack data at the cold extreme to estimate two parameters: TL, the lower temperature at which enzyme activity is 50% of the maximum value, and the Arrhenius constant associated with development below this threshold. In order to make a best guess of these parameters we write a very simple sum of squares function and find parameters that minimize it locally. For TL, we use an initial value that is a linear approximation of the 50% threshold from the data.

## Supplemental information: Solving the population model DDEs (Python code)

The code in DDE_one_pop.py solves the system of equations for the population model for a given set of parameter values. 
The code relies on PyDDE and numpy modules. 

### Parameter definitions

The function define_params defines all demographic parameters for three species, given as the only input to the function: spp which can be either 'trop', 'temp', or 'med'. The function returns an array of all parameter values.

There are constant parameter values that never vary in our scenarios: a year has 365 days, climate change occurs over 100 years, the time step of the system of equations is one day, and one of the dde solver parameter values that affects behaviour around zero has a specific values (note this has little effect).

In each climate change scenario we must define the species affected (spp), the time over which to solve (max_years), the number of years to plot or save (keep_years), the change in mean temperature (delta_mean), the change in the amplitude of seasonal fluctuations (delta_ampl), the functional form of the response of competition strength to temperature (q_form), the process that is density dependent (dd), and solver parameters (tol, dde_dt, dde_hbsize).

### Functions that define the thermal repsonses of life history traits

First we define the seasonal fluctuations in temperature as T(t). We then define the general functional forms for the Boltzmann-Arrhenius function, the Sharpe-Schoolfield function, and a general Gaussian relationship. We then define specific relationships for the thermal response of reproduction, mortality, development, and competition.

### Population model

The main function in the code ddegrad describes the system of delay differential equations that control the way the population state variables change through time. It takes an array of current state variables (s = {J, A, S_J, tau}), an empty array of parameters c, and the current time t. It returns an array {dJ/dt, dA/dt, dSJ/dt, dtau/dt}. The function also contains a clause to capture warnings for when the algorithm fails to accurately solve the system when one of the state variables approaches zero. If any of the state variables in s crosses to fall below zero (as should be impossible with non-negatve initial values), the function prints the current system variables. Becuase it prints at every dt at every time step, the solver slows down appreciably, alerting the user to the issue. To solve the misbehaviour around zero, the user may decrease tol and dde_dt. This may nessecitate increasing the size of the history buffer (hb_size) that can be accessed by the algorithm. The results are written to a txt file that can be read and plotted with your data visualization program of choice.



