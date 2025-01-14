dt = 1/24;
params.Arrival=1;
params.betamax=0.6;
params.N=20000;
params.I0=10;
params.VC=10;
params.IntT = dt;
params.xi = 2;
params.EndT = 250;
params.lam=0;
params.k=0;
params.ModelType = 1;
params.Ph=0.8;
params.Pl=0.2;
params.minP = 0.5;
params.changeP = 0.3;
params.shape = 2;
params.NoTest = 100;
params.BetaStar = 2/3*params.betamax;
params.FoC = 0.01; % Flu is 0.01, Covid is 0.1
params.IsolationDurtion = 7;

% TIV params
paramsT.Beta_I= [2.7*10^-5 3.2*10^-5];
paramsT.delta =[4 5.2];
paramsT.p = [1.2*10^-2 4.6*10^-2];
paramsT.c = [3 5.2];
paramsT.IntT = dt;
paramsT.T0 = 4*10^8;
paramsT.L0=0;
paramsT.I0=0;
paramsT.V0=[9.3*10^-2 7.5*10^-2];
paramsT.maxT = 15;
paramsT.k = 4;
paramsT.modelNum = 1;
params.VL_50 = 1000;

beta_tau = {@(x) params.betamax*x.^params.xi./(x.^params.xi+params.VL_50^params.xi)};

condition = @(x) x(2)==0;


model.paramsT = paramsT;
model.params=params;
model.condition=condition;
model.beta_tau=beta_tau;


% states_MC_SF_9 = {};
% states_MC_SF_11 = {};
% states_MC_SF_13 = {};
% states_MC_SF_15 = {};
% states_MC_SF_17 = {};
% 
% agent_TestTime_SF_9= {};
% agent_InfectRecoverTime_SF_9 = {};
% agent_TestTime_SF_11 = {};
% agent_InfectRecoverTime_SF_11 = {};
% agent_TestTime_SF_13 = {};
% agent_InfectRecoverTime_SF_13 = {};
% agent_TestTime_SF_15 = {};
% agent_InfectRecoverTime_SF_15 = {};
% agent_TestTime_SF_17 = {};
% agent_InfectRecoverTime_SF_17 = {};

Vstar = params.betamax*2/3;

for i = 1:100
    % 
    % [Ytime_H, states_H,InfectIdtau_H,agent_H,~,Isolation_H,DailyIsolation_H] = ABM_TIV_Policy(model,5,params.EndT,1,0,i);
    % 
    % states_MC_H{i}=states_H;
    % 
    % InfectList = max(find([agent_H(:).R]==1));
    % for I = 1:InfectList
    %     agent_TestTime_H{I,i} = agent_H(I).TestTime;
    %     agent_InfectRecoverTime_H{I,i}= [agent_H(I).InfectT agent_H(I).RecoverT];
    % end
    % 
    % IsolationList_H{i} = Isolation_H;
    % DailyIsolationList_H(i,:)=DailyIsolation_H;
    % CummulativeIso_H(i) = sum([agent_H.Isolation]);
    % 
    [Ytime_B, states_B,InfectIdtau_B,agent_B,~,Isolation_B,DailyIsolation_B] = ABM_TIV_Policy(model,5,params.EndT,1,1,i);

    states_MC_B{i}=states_B;

    % InfectList = max(find([agent_B(:).R]==1));
    % for I = 1:InfectList
    %     agent_TestTime_B{I,i} = agent_B(I).TestTime;
    %     agent_InfectRecoverTime_B{I,i}= [agent_B(I).InfectT agent_B(I).RecoverT];
    % end
    % 
    % IsolationList_B{i} = Isolation_B;
    % DailyIsolationList_B(i,:)=DailyIsolation_B;
    % CummulativeIso_B(i) = sum([agent_B.Isolation]);

    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    


    params.EndT = max(5,ceil(Peaktime/24)-35);
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,Isolation_TN_BL,DailyIsolation_BL] = ABM_TIV_Policy(model,5,params.EndT,1,1,i);
    Beta_inc__35{i}={};
    Beta_dec__35{i}={};
    InfectN = InfectIdtau_BL{end};
    for ii =1:length(InfectN)
        if agent_BL(InfectN(ii)).betta(params.EndT*24+1)>agent_BL(InfectN(ii)).betta(params.EndT*24)
                Beta_inc__35{i}=[Beta_inc__35{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];              
        else
            Beta_dec__35{i}=[Beta_dec__35{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];
        end    
    end
    Beta_inc__35{i} = cell2mat(Beta_inc__35{i});
    Beta_dec__35{i} = cell2mat(Beta_dec__35{i});

    params.EndT = 250;
    model.params=params;
    [Ytime_SF__35, states_SF__35,InfectIdtau_SF__35,agent_SF__35,~,Isolation_SF__35,DailyIsolation_SF__35] = ABM_TIV_Policy_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,min(5,ceil(Peaktime/24)-35),params.EndT,1,0,i,Isolation_TN_BL,DailyIsolation_BL);
    states_MC_SF__35{i}=states_SF__35;

    IsolationList_SF__35{i} = Isolation_SF__35;
    DailyIsolationList_SF__35(i,:)=DailyIsolation_SF__35;
    CummulativeIso_SF__35(i) = sum(DailyIsolation_SF__35);
    
    
    InfectList = max(find([agent_SF__35(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF__35{I,i} = agent_SF__35(I).TestTime;
        agent_InfectRecoverTime_SF__35{I,i}= [agent_SF__35(I).InfectT agent_SF__35(I).RecoverT];
    end

    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = max(5,ceil(Peaktime/24)-21);
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,Isolation_TN_BL,DailyIsolation_BL] = ABM_TIV_Policy(model,5,params.EndT,1,1,i);
    Beta_inc__21{i}={};
    Beta_dec__21{i}={};
    InfectN = InfectIdtau_BL{end};
    for ii =1:length(InfectN)
        if agent_BL(InfectN(ii)).betta(params.EndT*24+1)>agent_BL(InfectN(ii)).betta(params.EndT*24)
                Beta_inc__21{i}=[Beta_inc__21{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];              
        else
            Beta_dec__21{i}=[Beta_dec__21{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];
        end    
    end
    Beta_inc__21{i} = cell2mat(Beta_inc__21{i});
    Beta_dec__21{i} = cell2mat(Beta_dec__21{i});

    params.EndT = 250;
    model.params=params;
    [Ytime_SF__21, states_SF__21,InfectIdtau_SF__21,agent_SF__21,~,Isolation_SF__21,DailyIsolation_SF__21] = ABM_TIV_Policy_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)-21,params.EndT,1,0,i,Isolation_TN_BL,DailyIsolation_BL);
    states_MC_SF__21{i}=states_SF__21;
    IsolationList_SF__21{i} = Isolation_SF__21;
    DailyIsolationList_SF__21(i,:)=DailyIsolation_SF__21;
    CummulativeIso_SF__21(i) = sum(DailyIsolation_SF__21);

    InfectList = max(find([agent_SF__21(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF__21{I,i} = agent_SF__21(I).TestTime;
        agent_InfectRecoverTime_SF__21{I,i}= [agent_SF__21(I).InfectT agent_SF__21(I).RecoverT];
    end

    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24)-7;
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,Isolation_TN_BL,DailyIsolation_BL] = ABM_TIV_Policy(model,5,params.EndT,1,1,i);
    Beta_inc__7{i}={};
    Beta_dec__7{i}={};
    InfectN = InfectIdtau_BL{end};
    for ii =1:length(InfectN)
        if agent_BL(InfectN(ii)).betta(params.EndT*24+1)>agent_BL(InfectN(ii)).betta(params.EndT*24)
                Beta_inc__7{i}=[Beta_inc__7{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];              
        else
            Beta_dec__7{i}=[Beta_dec__7{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];
        end    
    end
    Beta_inc__7{i} = cell2mat(Beta_inc__7{i});
    Beta_dec__7{i} = cell2mat(Beta_dec__7{i});

    params.EndT = 250;
    model.params=params;
    [Ytime_SF__7, states_SF__7,InfectIdtau_SF__7,agent_SF__7,~,Isolation_SF__7,DailyIsolation_SF__7] = ABM_TIV_Policy_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)-7,params.EndT,1,0,i,Isolation_TN_BL,DailyIsolation_BL);
    states_MC_SF__7{i}=states_SF__7;
    IsolationList_SF__7{i} = Isolation_SF__7;
    DailyIsolationList_SF__7(i,:)=DailyIsolation_SF__7;
    CummulativeIso_SF__7(i) = sum(DailyIsolation_SF__7);
    
    InfectList = max(find([agent_SF__7(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF__7{I,i} = agent_SF__7(I).TestTime;
        agent_InfectRecoverTime_SF__7{I,i}= [agent_SF__7(I).InfectT agent_SF__7(I).RecoverT];
    end

    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24)+7;
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,Isolation_TN_BL,DailyIsolation_BL] = ABM_TIV_Policy(model,5,params.EndT,1,1,i);
    Beta_inc_7{i}={};
    Beta_dec_7{i}={};
    InfectN = InfectIdtau_BL{end};
    for ii =1:length(InfectN)
        if agent_BL(InfectN(ii)).betta(params.EndT*24+1)>agent_BL(InfectN(ii)).betta(params.EndT*24)
                Beta_inc_7{i}=[Beta_inc_7{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];              
        else
            Beta_dec_7{i}=[Beta_dec_7{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];
        end    
    end
    Beta_inc_7{i} = cell2mat(Beta_inc_7{i});
    Beta_dec_7{i} = cell2mat(Beta_dec_7{i});

    params.EndT = 250;
    model.params=params;
    [Ytime_SF_7, states_SF_7,InfectIdtau_SF_7,agent_SF_7,~,Isolation_SF_7,DailyIsolation_SF_7] = ABM_TIV_Policy_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)+7,params.EndT,1,0,i,Isolation_TN_BL,DailyIsolation_BL);
    states_MC_SF_7{i}=states_SF_7;

    IsolationList_SF_7{i} = Isolation_SF_7;
    DailyIsolationList_SF_7(i,:)=DailyIsolation_SF_7;
    CummulativeIso_SF_7(i) = sum(DailyIsolation_SF_7);
    
    InfectList = max(find([agent_SF_7(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF_7{I,i} = agent_SF_7(I).TestTime;
        agent_InfectRecoverTime_SF_7{I,i}= [agent_SF_7(I).InfectT agent_SF_7(I).RecoverT];
    end

    

    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24)+21;
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,Isolation_TN_BL,DailyIsolation_BL] = ABM_TIV_Policy(model,5,params.EndT,1,1,i);
    Beta_inc_21{i}={};
    Beta_dec_21{i}={};
    InfectN = InfectIdtau_BL{end};
    for ii =1:length(InfectN)
        if agent_BL(InfectN(ii)).betta(params.EndT*24+1)>agent_BL(InfectN(ii)).betta(params.EndT*24)
                Beta_inc_21{i}=[Beta_inc_21{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];              
        else
            Beta_dec_21{i}=[Beta_dec_21{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];
        end    
    end
    Beta_inc_21{i} = cell2mat(Beta_inc_21{i});
    Beta_dec_21{i} = cell2mat(Beta_dec_21{i});
    params.EndT = 250;
    model.params=params;
    [Ytime_SF_21, states_SF_21,InfectIdtau_SF_21,agent_SF_21,~,Isolation_SF_21,DailyIsolation_SF_21] = ABM_TIV_Policy_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)+21,params.EndT,1,0,i,Isolation_TN_BL,DailyIsolation_BL);
    states_MC_SF_21{i}=states_SF_21;

    IsolationList_SF_21{i} = Isolation_SF_21;
    DailyIsolationList_SF_21(i,:)=DailyIsolation_SF_21;
    CummulativeIso_SF_21(i) = sum(DailyIsolation_SF_21);
    
    InfectList = max(find([agent_SF_21(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF_21{I,i} = agent_SF_21(I).TestTime;
        agent_InfectRecoverTime_SF_21{I,i}= [agent_SF_21(I).InfectT agent_SF_21(I).RecoverT];
    end

    Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    params.EndT = ceil(Peaktime/24)+35;
    model.params=params;
    [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,Isolation_TN_BL,DailyIsolation_BL] = ABM_TIV_Policy(model,5,params.EndT,1,1,i);
    Beta_inc_35{i}={};
    Beta_dec_35{i}={};
    InfectN = InfectIdtau_BL{end};
    for ii =1:length(InfectN)
        if agent_BL(InfectN(ii)).betta(params.EndT*24+1)>agent_BL(InfectN(ii)).betta(params.EndT*24)
                Beta_inc_35{i}=[Beta_inc_35{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];              
        else
            Beta_dec_35{i}=[Beta_dec_35{i}, agent_BL(InfectN(ii)).betta(params.EndT*24+1)];
        end    
    end
    Beta_inc_35{i} = cell2mat(Beta_inc_35{i});
    Beta_dec_35{i} = cell2mat(Beta_dec_35{i});
    params.EndT = 250;
    model.params=params;
    [Ytime_SF_35, states_SF_35,InfectIdtau_SF_35,agent_SF_35,~,Isolation_SF_35,DailyIsolation_SF_35] = ABM_TIV_Policy_Conti(model,...
        Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)+35,params.EndT,1,0,i,Isolation_TN_BL,DailyIsolation_BL);
    states_MC_SF_35{i}=states_SF_35;

    IsolationList_SF_35{i} = Isolation_SF_35;
    DailyIsolationList_SF_35(i,:)=DailyIsolation_SF_35;
    CummulativeIso_SF_35(i) = sum(DailyIsolation_SF_35);

    InfectList = max(find([agent_SF_35(:).R]==1));
    for I = 1:InfectList
        agent_TestTime_SF_35{I,i} = agent_SF_35(I).TestTime;
        agent_InfectRecoverTime_SF_35{I,i}= [agent_SF_35(I).InfectT agent_SF_35(I).RecoverT];
    end
    % Peaktime = max(find(states_B(:,2)==max(states_B(:,2))));
    % params.EndT = ceil(Peaktime/24)+17;
    % model.params=params;
    % [Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL] = ABM_TIV_Policy(model,5,params.EndT,1,1,i);
    % 
    % params.EndT = 200;
    % model.params=params;
    % [Ytime_SF_17, states_SF_17,InfectIdtau_SF_17,agent_SF_17] = ABM_TIV_Policy_Conti(model,...
    %     Ytime_BL, states_BL,InfectIdtau_BL,agent_BL,Isolation_BL,ceil(Peaktime/24)+17,params.EndT,1,0,i);
    % states_MC_SF_17{i}=states_SF_17;
    % 
    % InfectList = max(find([agent_SF_17(:).R]==1));
    % for I = 1:InfectList
    %     agent_TestTime_SF_17{I,i} = agent_SF_17(I).TestTime;
    %     agent_InfectRecoverTime_SF_17{I,i}= [agent_SF_17(I).InfectT agent_SF_17(I).RecoverT];
    % end

    clear agent_H agent_B agent_SF__35 agent_SF__21 agent_SF__7 agent_SF_7 agent_SF_21 agent_SF_35

end

% save("MC100_MT_DifferetTime_Final_v2.mat",...
%     'states_MC_B','agent_TestTime_B','agent_InfectRecoverTime_B',...
%     'states_MC_SF__21','agent_TestTime_SF__21','agent_InfectRecoverTime_SF__21',...
%     'states_MC_SF__35','agent_TestTime_SF__35','agent_InfectRecoverTime_SF__35',...
%     'states_MC_SF__7','agent_TestTime_SF__7','agent_InfectRecoverTime_SF__7',...
%     'states_MC_B','agent_TestTime_B','agent_InfectRecoverTime_B',...
%     'states_MC_SF_7','agent_TestTime_SF_7','agent_InfectRecoverTime_SF_7',...
%     'states_MC_SF_35','agent_TestTime_SF_35','agent_InfectRecoverTime_SF_35',...
%     'states_MC_SF_21','agent_TestTime_SF_21','agent_InfectRecoverTime_SF_21')
% 
% save("MC100_BetaVlue.mat",'Beta_inc__35','Beta_dec__35','Beta_inc__21','Beta_dec__21',...
%     'Beta_inc__7','Beta_dec__7','Beta_inc_35','Beta_dec_35','Beta_inc_21','Beta_dec_21',...
%     'Beta_inc_7','Beta_dec_7')

save("MC100_MT_isolation_DiffTime.mat",'states_MC_SF__21','states_MC_SF__35','states_MC_SF__7',...
    'states_MC_SF_7','states_MC_SF_35','states_MC_SF_21',...
    "IsolationList_SF__35","DailyIsolationList_SF__35","CummulativeIso_SF__35",...
    "IsolationList_SF__21","DailyIsolationList_SF__21","CummulativeIso_SF__21",...
    "IsolationList_SF__7","DailyIsolationList_SF__7","CummulativeIso_SF__7",...
    "IsolationList_SF_7","DailyIsolationList_SF_7","CummulativeIso_SF_7",...
    "IsolationList_SF_21","DailyIsolationList_SF_21","CummulativeIso_SF_21",...
    "IsolationList_SF_35","DailyIsolationList_SF_35","CummulativeIso_SF_35")