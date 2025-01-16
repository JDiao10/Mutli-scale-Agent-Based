# Multi-scale agent-based model and policy-making
## The Basic_code folder is the SIR-TIV agent-based model, including TIV.m and AMB_TIV.m
TIV.m is the basic target-limited cell (TIV) model that can generate the state of T, I, V at different time steps;

ABM_TIV.m is the agent-based model that needs to have the parameters to run;

## Policy_making
To implement different testing-isolation strategies (isolate all, isolate high, adaptive) and save the simulated data, please read the readme in Policy_making for more details.

## Plot_Code
Generating figures based on the simulated data from Policy_making.
