# Discrete Time Heterogenous Agent Model

The algorithms used to solve the model were in large part
originally written by Greg Kaplan and modified by Brian Livingston
(livingstonb@uchicago.edu).
The *aux_lib* directory contains code provided by Mario Miranda and Paul Fackler
via the CompEcon toolbox. In several other cases, this repository
contains code written by others, with citations where possible.

## Using the master script

The code is ultimately executed from the script *master.m*.
From this script, first set options by assigning values to the *runopts* structure in the OPTIONS section.
Note that the all of the computational code is called in *main.m*.
The master script just serves as a wrapper for that file.

### Parameters

Within the default parameters script, a parameterization is assigned to each structure within a structure array. The easiest way to select a specific parameterization to run is to set the *number* field of runopts equal to the index of the desired parameterization within the structure array.
To run your desired parameterization, either modify *code/+setup/parameters.m* or create a new parameters file.
All parameter defaults are set in the class file *code/+setup/Params.m* and any values set in the selected parameters file will override the default value of the given variable.
The *make_adjustments* method of the Params class automatically adjusts some parameter values based on other parameter values so be careful not to make those adjustments manually,
e.g. if the model frequency is set to quarterly then the discount rate (which is assumed to be annualized)
is converted to a quarterly rate.
See the *make_adjustments* methods for details.
If the user needs to wrap the model in a non-linear solver to match user-specified moments, this can be done in the parameters file. See the next section for details.

### Calibration

The class defined in *code/+solver/DHACalibrator.m* is used to assist with matching moments. This class is implemented by creating a DHACalibrator instance in the parameters file, which is then assigned to an attribute of the Params object returned by the parameters file. This object allows the user to easily set the parameters that
need to be calibrated as well as the moments to match. See the class file for details.
Note: if convergence fails, betaH0 and/or betaL may need to be adjusted.
betaL is the lower bound picked for beta during iteration and
betaH0 is the adjustment factor to the upper bound. The code will
guess a theoretical upper bound, and then will add betaH0 to
to that value.

### Output

The output table can is automatically produced upon completion of the code.
Results are stored in the 'results' structure. Its 'direct' property
contains results found from computing the stationary distribution
using non-simulation numerical methods. The 'sim' property contains results
found from simulation, if the option is turned on.

## Running code on the server (Midway)

### Table creation