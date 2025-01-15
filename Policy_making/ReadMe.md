# Agent-based model that can implement the testing-isolation strategies
## The implementation code
ABM_TIV_Policy.m is the general code for applying different strategies, in which the testing-isolation strategy starts at the duration based on the input and until the ending time.

ABM_TIV_Policy_Conti.m is the code that can based on the output from ABM_TIV_Policy.m and then change different strategies until the end time of the outbreak. (Apply for changing the policy from isolate-all to adaptive strategy)

ABM_TIV_Policy_sym.m and ABM_TIV_Policy_sym_Conti.m are the same operation as mentioned above but include some infected agents that are asymptomatic.

## Codes produce the simulated data.
