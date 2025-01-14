dt = 1/24;
params.Arrival=1;
params.betamax=0.3;
params.N=20000;
params.I0=10;
params.VC=10;
params.IntT = dt;
params.xi = 2;
params.EndT = 250;
params.lam=0;
params.k=0;
params.ModelType = 2;
params.minP = 0.5;
params.changeP = 0.3;
params.NoTest = 100;
params.BetaStar=params.betamax*2/3;
params.Psym = 0.7;
params.Recover = 10^-1;

% TIV params
paramsT.Beta_I= 4.71*10^-8;
paramsT.delta =1.07;
paramsT.p = 3.07;
paramsT.c = 2.40;
paramsT.IntT = dt;
paramsT.T0 = 4*10^8;
paramsT.L0=0;
paramsT.I0=0;
paramsT.maxT = 21;
paramsT.k = 4;
paramsT.modelNum = 1;
params.VL_50 = 1000;

params.TimeTest = randi([1 24],params.N,paramsT.maxT);

beta_tau = {@(x) params.betamax*x^params.xi/(x^params.xi+params.VL_50^params.xi)};

condition = @(x) x(2)==0;

model.paramsT = paramsT;
model.params=params;
model.condition=condition;
model.beta_tau=beta_tau;



states_MC ={};
agent_TestTime = {};
agent_InfectRecoverTime = {};
timeX_MC = {};
Beta_Full = {};
Beta_Begin = {};
Beta_Recover = {};
Beta_BeginB={};
Beta_BeginS={};
Beta_RecoverB={};
Beta_RecoverS={};
agent_sym_H = {};
agent_sym_B ={};
agent_sym_SF ={};
Vstar = params.betamax*2/3;
% isolate all first and then isolate high only


for i = 1:20
    params.EndT = 250;
    model.params = params;
    [Ytime_H, states_H,InfectIdtau_H,agent_H,~] = ABM_TIV_Policy_sym(model,5,params.EndT,1,0,i);
    
    states_MC_H{i}=states_H;

    StepTime = 1/24;
    
    InfectList = max(find([agent_H(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_H{I,i} = agent_H(I).TestTime;
        agent_InfectRecoverTime_H{I,i}= [agent_H(I).InfectT agent_H(I).RecoverT];
        agent_sym_H{I,i} = agent_H(I).symtoms;
    end
    
    
    params.EndT= 250;
    model.params=params;
    [Ytime_B, states_B,InfectIdtau_B,agent_B,~] = ABM_TIV_Policy_sym(model,5,params.EndT,1,1,i);
    
    states_MC_B{i}=states_B;
        
        InfectList = max(find([agent_B(:).R]==1));
        for I = 1:InfectList
            agent_TestTime_B{I,i} = agent_B(I).TestTime;
            agent_InfectRecoverTime_B{I,i}= [agent_B(I).InfectT agent_B(I).RecoverT];
            agent_sym_B{I,i} = agent_B(I).symtoms;
        end
    
    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24);
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL] = ABM_TIV_Policy_sym(model,5,params.EndT,1,1,i);
    
    params.EndT= 250;
    model.params=params;
    [Ytime_SF, states_SF,InfectIdtau_SF,agent_SF] = ABM_TIV_Policy_sym_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24),params.EndT,1,0,i);

    states_MC_SF{i}=states_SF;
        
    InfectList = max(find([agent_SF(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF{I,i} = agent_SF(I).TestTime;
        agent_InfectRecoverTime_SF{I,i}= [agent_SF(I).InfectT agent_SF(I).RecoverT];
        agent_sym_SF{I,i} = agent_SF(I).symtoms;
    end

    clear agent_H agent_B agent_SF
end

save("MC_Covid_Data_PolicyH_MT_sym.mat","states_MC_H","agent_InfectRecoverTime_H","agent_TestTime_H",'agent_sym_H')

save("MC_Covid_Data_PolicyB_MT_sym.mat",'states_MC_B','agent_TestTime_B','agent_InfectRecoverTime_B','agent_sym_B')

save("Policydata_Covid_StopL_MT.mat_sym",'states_MC_SF','agent_TestTime_SF','agent_InfectRecoverTime_SF','agent_sym_SF')




