clear all
dt = 1/24;
params.Arrival=1;
params.betamax=0.6;
params.N=20000;
params.I0=10;
params.T=60;
params.VC=10;
params.IntT = dt;
params.xi = 2;

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

beta_tau = {@(x) params.betamax*x^params.xi/(x^params.xi+params.VL_50^params.xi)};

condition = @(x) x(2)==0;


model.paramsT = paramsT;
model.params=params;
model.condition=condition;
model.beta_tau=beta_tau;
tic
[Ytime,States,InfectIdtau,agent]=ABM_TIV(model);
toc

InfectList = max(find([agent(:).R]==1));
beta_Begin = [];
beta_Recover = [];
StepTime = 1/24;
for i = 1:InfectList
    MaxbetaTime = max(find(agent(i).betta==max(agent(i).betta)));
    beta_Begin = [beta_Begin, [agent(i).InfectT:StepTime:MaxbetaTime*StepTime-StepTime;...
        agent(i).betta(single(agent(i).InfectT/StepTime+1):MaxbetaTime)]];
    beta_Recover = [beta_Recover, [MaxbetaTime*StepTime+StepTime:singleStepTime:agent(i).RecoverT;...
        agent(i).betta(MaxbetaTime+1:single(agent(i).RecoverT/StepTime))]];
end
beta_Begin = sortrows(beta_Begin',[1])';
beta_Recover = sortrows(beta_Recover',[1])';
timeX=0:StepTime:max(beta_Recover(1,:));

save("ABM_TIV_value0_6_20000t1hr_new_DiffTIV.mat",'agent','timeX','beta_Begin','beta_Recover','States')
