dt = 1/24;
params.Arrival=1;
params.betamax=0.6;
params.N=20000;
params.I0=10;
params.VC=10;
params.IntT = dt;
params.xi = 2;
params.EndT = 220;
params.lam=0;
params.k=0;
params.ModelType = 2;
params.minP = 0.5;
params.changeP = 0.3;
params.BetaStar = params.betamax*2/3;
params.Psym = 0.6;
params.Recover = 10^-2;

% TIV params
paramsT.Beta_I= [2.7*10^-5 3.2*10^-5];
paramsT.delta =[4 5.2];
paramsT.p = [1.2*10^-2 4.6*10^-2];
paramsT.c = [3 5.2];
paramsT.IntT = dt;
paramsT.T0 = 4*10^8;
paramsT.L0=0;
paramsT.I0=0;
paramsT.maxT = 15;
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


states_MC_H ={};
states_MC_B={};
states_MC_SF = {};
states_MC_SF_1 = {};
states_MC_SF_3 = {};
states_MC_SF_5 = {};
states_MC_SF_7 = {};
agent_TestTime_H = {};
agent_InfectRecoverTime_H = {};
agent_TestTime_B = {};
agent_InfectRecoverTime_B = {};
agent_TestTime_SF = {};
agent_InfectRecoverTime_SF = {};
agent_TestTime_SF_1= {};
agent_InfectRecoverTime_SF_1 = {};
agent_TestTime_SF_3 = {};
agent_InfectRecoverTime_SF_3 = {};
agent_TestTime_SF_5 = {};
agent_InfectRecoverTime_SF_5 = {};
agent_TestTime_SF_7 = {};
agent_InfectRecoverTime_SF_7 = {};
agent_sym_H = {};
agent_sym_B ={};
agent_sym_SF ={};
agent_sym_SF_1 ={};
agent_sym_SF_3 ={};
agent_sym_SF_5 ={};
agent_sym_SF_7 ={};
Vstar = params.betamax*2/3;


for i = 1:20

    [Ytime_H, states_H,InfectIdtau_H,agent_H,~] = ABM_TIV_Policy_sym(model,5,params.EndT,1,0,i);
    
    
    states_MC_H{i}=states_H;
    
     StepTime = 1/24;
        
        InfectList = max(find([agent_H(:).R]==1));
        for I = 1:InfectList
            agent_TestTime_H{I,i} = agent_H(I).TestTime;
            agent_InfectRecoverTime_H{I,i}= [agent_H(I).InfectT agent_H(I).RecoverT];
            agent_sym_H{I,i} = agent_H(I).symtoms;
        end
            
    
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
    
    params.EndT = 220;
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

    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24)+1;
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL] = ABM_TIV_Policy_sym(model,5,params.EndT,1,1,i);
    
    params.EndT = 220;
    model.params=params;
    [Ytime_SF_1, states_SF_1,InfectIdtau_SF_1,agent_SF_1] = ABM_TIV_Policy_sym_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)+1,params.EndT,1,0,i);
    states_MC_SF_1{i}=states_SF_1;
    
    InfectList = max(find([agent_SF_1(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF_1{I,i} = agent_SF_1(I).TestTime;
        agent_InfectRecoverTime_SF_1{I,i}= [agent_SF_1(I).InfectT agent_SF_1(I).RecoverT];
        agent_sym_SF_1{I,i} = agent_SF_1(I).symtoms;
    end
    
    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24)+3;
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL] = ABM_TIV_Policy_sym(model,5,params.EndT,1,1,i);
    
    params.EndT = 220;
    model.params=params;
    [Ytime_SF_3, states_SF_3,InfectIdtau_SF_3,agent_SF_3] = ABM_TIV_Policy_sym_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)+3,params.EndT,1,0,i);
    states_MC_SF_3{i}=states_SF_3;
    
    InfectList = max(find([agent_SF_3(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF_3{I,i} = agent_SF_3(I).TestTime;
        agent_InfectRecoverTime_SF_3{I,i}= [agent_SF_3(I).InfectT agent_SF_3(I).RecoverT];
        agent_sym_SF_3{I,i} = agent_SF_3(I).symtoms;
    end
    
    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24)+5;
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL] = ABM_TIV_Policy_sym(model,5,params.EndT,1,1,i);
    
    params.EndT = 220;
    model.params=params;
    [Ytime_SF_5, states_SF_5,InfectIdtau_SF_5,agent_SF_5] = ABM_TIV_Policy_sym_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)+5,params.EndT,1,0,i);
    states_MC_SF_5{i}=states_SF_5;
    
    InfectList = max(find([agent_SF_5(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF_5{I,i} = agent_SF_5(I).TestTime;
        agent_InfectRecoverTime_SF_5{I,i}= [agent_SF_5(I).InfectT agent_SF_5(I).RecoverT];
        agent_sym_SF_5{I,i} = agent_SF_5(I).symtoms;
    end
    
    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24)+7;
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL] = ABM_TIV_Policy_sym(model,5,params.EndT,1,1,i);
    
    params.EndT = 220;
    model.params=params;
    [Ytime_SF_7, states_SF_7,InfectIdtau_SF_7,agent_SF_7] = ABM_TIV_Policy_sym_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)+7,params.EndT,1,0,i);
    states_MC_SF_7{i}=states_SF_7;
    
    InfectList = max(find([agent_SF_7(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF_7{I,i} = agent_SF_7(I).TestTime;
        agent_InfectRecoverTime_SF_7{I,i}= [agent_SF_7(I).InfectT agent_SF_7(I).RecoverT];
        agent_sym_SF_7{I,i} = agent_SF_7(I).symtoms;
    end
    
    
    clear agent_H agent_B agent_SF agent_SF_1 agent_SF_3 agent_SF_5 agent_SF_7

end

save("MC_Flu_Data_PolicyH_MT_sym.mat","states_MC_H","agent_InfectRecoverTime_H","agent_TestTime_H",'agent_sym_H')

save("MC_Flu_Data_PolicyB_MT_sym.mat",'states_MC_B','agent_TestTime_B','agent_InfectRecoverTime_B','agent_sym_B')

save("Policydata_Flu_StopL_MT_sym.mat",'states_MC_SF','agent_TestTime_SF','agent_InfectRecoverTime_SF','agent_sym_SF')

save("MC_MT_DifferetTime_sym.mat",'states_MC_SF','agent_TestTime_SF','agent_InfectRecoverTime_SF','agent_sym_SF',...
    'states_MC_SF_1_sym','agent_TestTime_SF_1','agent_InfectRecoverTime_SF_1','agent_sym_SF_1',...
    'states_MC_SF_3_sym','agent_TestTime_SF_3','agent_InfectRecoverTime_SF_3','agent_sym_SF_3',...
    'states_MC_SF_5_sym','agent_TestTime_SF_5','agent_InfectRecoverTime_SF_5','agent_sym_SF_5',...
    'states_MC_SF_7_sym','agent_TestTime_SF_7','agent_InfectRecoverTime_SF_7','agent_sym_SF_7')






