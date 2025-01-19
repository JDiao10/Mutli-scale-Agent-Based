# Code for generating results of the paper "Should public health policy exempt cases with low viral load from isolation during an epidemic?: a modelling study". 
## The Basic_code folder is the SIR-TIV agent-based model, including TIV.m and AMB_TIV.m
TIV.m is the basic target-limited cell (TIV) model that can generate the state of T, I, V at different time steps;

ABM_TIV.m is the agent-based model that includes the parameters of the TIV within-host model, the parameters at the population level, and more details in the readme in Basic_code;

## Policy_making
To implement different testing-isolation strategies (isolate all, isolate high, adaptive) and save the simulated data, please read the readme in Policy_making for more details.

## Plot_Code
Generating figures based on the simulated data from Policy_making.
