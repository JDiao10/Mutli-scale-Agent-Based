dt = 1/24;
params.Arrival=1;
params.betamax=0.3;
params.N=20000;
params.I0=10;
params.VC=10;
params.IntT = dt;
params.xi = 2;
params.EndT = 200;
params.lam=0;
params.k=0;
params.ModelType = 2;
params.Ph=0.5;
params.Pl=0.1;
params.minP = 0.5;
params.changeP = 0.3;
params.NoTest = 1;
params.BetaStar = 2/3*params.betamax;
params.FoC = 0.01; % Flu is 0.01, Covid is 0.1
params.IsolationDurtion = 14;

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

beta_tau = {@(x) params.betamax*x.^params.xi./(x.^params.xi+params.VL_50^params.xi)};


condition = @(x) x(2)==0;


model.paramsT = paramsT;
model.params=params;
model.condition=condition;
model.beta_tau=beta_tau;

states_MC_No={};
states_MC_H ={};
states_MC_B={};
states_MC_SF = {};
agent_TestTime_No = {};
agent_InfectRecoverTime_No= {};
agent_InfectRecoverTime_H = {};
agent_TestTime_B = {};
agent_InfectRecoverTime_B = {};
agent_TestTime_SF = {};
agent_InfectRecoverTime_SF = {};
agent_sym_H = {};
agent_sym_B ={};
agent_sym_SF ={};
Vstar = params.betamax*2/3;


for i = 1:100
    % [Ytime_No, states_No,InfectIdtau_No,agent_No,~] = ABM_TIV_Policy(model,0,0,0,0,i);
    % 
    % states_MC_No{i}=states_No;
    % 
    % StepTime = 1/24;
    % 
    % InfectList = max(find([agent_No(:).R]==1));
    % for I = 1:InfectList
    %     agent_TestTime_No{I,i} = agent_No(I).TestTime;
    %     agent_InfectRecoverTime_No{I,i}= [agent_No(I).InfectT agent_No(I).RecoverT];
    % end
    % 
    % 
    % 
    [Ytime_H, states_H,InfectIdtau_H,agent_H,~,Isolation_H,DailyIsolation_H] = ABM_TIV_Policy(model,10,params.EndT,1,0,i);

    states_MC_H{i}=states_H;

    StepTime = 1/24;

    InfectList = max(find([agent_H(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_H{I,i} = agent_H(I).TestTime;
        agent_InfectRecoverTime_H{I,i}= [agent_H(I).InfectT agent_H(I).RecoverT];
    end

    IsolationList_H{i} = Isolation_H;
    DailyIsolationList_H(i,:)=DailyIsolation_H;
    CummulativeIso_H(i) = sum(DailyIsolation_H);  

    [Ytime_B, states_B,InfectIdtau_B,agent_B,~,Isolation_B,DailyIsolation_B] = ABM_TIV_Policy(model,10,params.EndT,1,1,i);

    states_MC_B{i}=states_B;
    
    InfectList = max(find([agent_B(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_B{I,i} = agent_B(I).TestTime;
        agent_InfectRecoverTime_B{I,i}= [agent_B(I).InfectT agent_B(I).RecoverT];
        agent_SIbeta{I,i} = agent_B(I).betta;
    end

    IsolationList_B{i} = Isolation_B;
    DailyIsolationList_B(i,:)=DailyIsolation_B;
    CummulativeIso_B(i) = sum(DailyIsolation_B);

    Beta_incST{i}={};
    Beta_decST{i}={};
    NoIST(i) = states_MC_B{i}(end,3);
    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    peakdayST(i) = ceil(max(find(states_MC_B{i}(:,2)==max(states_MC_B{i}(:,2))))/24);

    for j = 1:NoIST(i)
        if agent_InfectRecoverTime_B{j,i}(1)<=peakdayST(i) & agent_InfectRecoverTime_B{j,i}(2)>=peakdayST(i)
            if agent_SIbeta{j,i}(peakdayST(i)*24+1)>agent_SIbeta{j,i}(peakdayST(i)*24)
                Beta_incST{i}=[Beta_incST{i} agent_SIbeta{j,i}(peakdayST(i)*24+1)];
            else
                Beta_decST{i}=[Beta_decST{i} agent_SIbeta{j,i}(peakdayST(i)*24+1)];
            end

        end
    end


    params.EndT = ceil(Peaktime/24);
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,Isolation_TN_BL,DailyIsolation_BL] = ABM_TIV_Policy(model,10,params.EndT,1,1,i);
    
    params.EndT = 200;
    model.params=params;
    [Ytime_SF, states_SF,InfectIdtau_SF,agent_SF,~,Isolation_SF,DailyIsolation_SF] = ABM_TIV_Policy_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24),params.EndT,1,0,i,Isolation_TN_BL,DailyIsolation_BL);
    states_MC_SF{i}=states_SF;
    
    InfectList = max(find([agent_SF(:).R]==1));
    IsolationList_SF{i} = Isolation_SF;
    DailyIsolationList_SF(i,:)=DailyIsolation_SF;
    CummulativeIso_SF(i) = sum(DailyIsolation_SF);
    for I = 1:InfectList
        agent_TestTime_SF{I,i} = agent_SF(I).TestTime;
        agent_InfectRecoverTime_SF{I,i}= [agent_SF(I).InfectT agent_SF(I).RecoverT];
    end
    Beta_incST{i} = cell2mat(Beta_incST{i});
    Beta_decST{i} = cell2mat(Beta_decST{i});
    clear clear agent_No agent_H agent_B agent_SF
end

%save("MC100_Covid_Data_No.mat","states_MC_No","agent_InfectRecoverTime_No","agent_TestTime_No")

%save("MC100_Covid_Data_PolicyH.mat","states_MC_H","agent_InfectRecoverTime_H","agent_TestTime_H")

%save("MC100_Covid_Data_PolicyB.mat",'states_MC_B','agent_TestTime_B','agent_InfectRecoverTime_B','Beta_incST','Beta_decST')

%save("MC100_Policydata_Covid_StopL.mat",'states_MC_SF','agent_TestTime_SF','agent_InfectRecoverTime_SF')

save("MC100_Covid_Data_Policy_Isolation_ST.mat","states_MC_H","states_MC_B","states_MC_SF",...
    "IsolationList_H","IsolationList_B","IsolationList_SF","DailyIsolationList_H",...
    "DailyIsolationList_B","DailyIsolationList_SF","CummulativeIso_H","CummulativeIso_B","CummulativeIso_SF")
