function [Ytime,states,InfectIdtau,agent,Isolation]=ABM_TIV_Policy(model,IntStartT,DurationT,OnlyH,HwithL,seed)
rng(seed)
beta_tau=model.beta_tau;
Get_beta=beta_tau{1};
condition=model.condition;
params = model.params;
I0 = params.I0;
TimeStep = params.IntT;
paramsT=model.paramsT;
tspan =0:paramsT.IntT:paramsT.maxT;
V0 = 10.^(2*rand(1,params.N));
Ph = params.Ph;
Pl = params.Pl;
minP = params.minP;
range = params.changeP;
ModelType = params.ModelType;
NoTest = params.NoTest;

% Record the person is infected at time \tau

for i = 1:params.N
    agent(i).I = 0;
    agent(i).S = 1;
    agent(i).R = 0;
    agent(i).InfectT=0;
    agent(i).RecoverT=0;
    agent(i).TIVs = {};
    agent(i).IID = [];
    agent(i).Test = NoTest;
    agent(i).Ptest = zeros(params.EndT/TimeStep+1,1);
    agent(i).Pposit = zeros(params.EndT/TimeStep+1,1);
    agent(i).TestTime = [];
end

InfectN=1:I0;

for i = 1:length(InfectN)
    agent(i).I=1;
    agent(i).S=0;
    [stateI] = TIV(model,V0(i));    
    %beta_Begin=[beta_Begin, [0;agent(i).betta(1)]];
    RecoverT = recover_fun(stateI,tspan,params.FoC);
    agent(i).betta = Get_beta(stateI(1:single(RecoverT/params.IntT)+1,3))';
    agent(i).RecoverT = RecoverT;   
    agent(i).TIVs = stateI;
    agent(i).IID = 0;
    agent(i).taulist = 0:TimeStep:RecoverT;
    taulist = agent(i).taulist;
    agent(i).Ptest(1:length(taulist)) =...
        Prob_get_test(agent(i).taulist,ModelType,Ph,Pl,TimeStep);
    agent(i).Pposit(1:length(taulist)) = ...
        Prob_testPosit_VL(agent(i).TIVs(1:length(taulist),3),minP,range);

end

InfectIdtau{1} = InfectN;

states = [sum([agent(:).S]) sum([agent(:).I]) sum([agent(:).R])];
% Record the Susceptiable individual number
SuscepN=I0+1:params.N;
% Record the time
Ytime=0;

Isolation = [];
%IntPre = params.IntPre;

InfectN = InfectIdtau{end};
SuscepN=max(InfectIdtau{end})+1:params.N;


% Record the Susceptiable individual number

% Record the time

        
  
for ii = 1:params.EndT/TimeStep
    tt=ii*TimeStep;

% Need to change this part (first need to get test
    if tt <= DurationT & tt>=IntStartT
        TestInfectN=InfectN(find([agent(InfectN).Test]>0));
        for i = 1:length(TestInfectN)
            % Infected hosts go for test
            if rand<agent(TestInfectN(i)).Ptest(ii)
                agent(TestInfectN(i)).Test = agent(TestInfectN(i)).Test-1;
                agent(TestInfectN(i)).TestTime = [agent(TestInfectN(i)).TestTime [tt;0]];
                % Test positive
                if rand<agent(TestInfectN(i)).Pposit(ii)
                    agent(TestInfectN(i)).TestTime(2,end)=1;
                    if isolate_policy(OnlyH,HwithL)==1
                        if agent(TestInfectN(i)).betta(ii)>=params.BetaStar
                            Isolation = [Isolation TestInfectN(i)];
                        end
                    elseif isolate_policy(OnlyH,HwithL) == 2  
                            Isolation = [Isolation TestInfectN(i)];                     
                    elseif isolate_policy(OnlyH,HwithL) == 3
                        if  agent(TestInfectN(i)).betta(ii)<params.BetaStar
                            Isolation = [Isolation TestInfectN(i)];
                        end
                    end
                end

            end
            
        end
        % agentBetta=[];
        % for i = 1:length(NoTestInfectN)
        %     agentBetta(i) = agent(NoTestInfectN(i)).betta(single(max(Ytime)/TimeStep)+1);
        % end
        % LowVl = find(agentBetta<=0.1);
        % if length(LowVl)>=1
        %     Isolation = [Isolation NoTestInfectN(LowVl(randperm(length(LowVl),round(IntPre*length(LowVl)))))];
        %     for j = 1:length(LowVl)
        %         agent(NoTestInfectN(LowVl(j))).Test= 0;
        %     end
        % end
    end


    
    % Isolated infected host to recover
    InfectN(ismember(InfectN,Isolation)) = [];
    if length(Isolation)>=1
        RecoverIso = [];
        for i = 1:length(Isolation)
            if tt>agent(Isolation(i)).RecoverT
                agent(Isolation(i)).I = 0;
                agent(Isolation(i)).R = 1;   
                RecoverIso = [RecoverIso i];
            end            
        end
        Isolation(RecoverIso)=[];
    end



    if ~condition(states(end,:))
        InfectPeople = []; % Record the infected people being recover for each time step
        for j = 1:length(InfectN)
            if tt<=agent(InfectN(j)).RecoverT
                
                
                if length(SuscepN)>0
                    ContactN = round(poissrnd(params.Arrival*TimeStep));
                    Contact = binornd(ContactN,length(SuscepN)/params.N);

                    if Contact>0
                        SuscepPeople=[];
                        for i = 1:Contact
                            
                            p_ill = agent(InfectN(j)).betta(single(tt/params.IntT));
        
                            p = rand;
        
                            if p<=p_ill                              
                                agent(SuscepN(i)).InfectT=tt;
                                agent(SuscepN(i)).I=1;                                
                                agent(SuscepN(i)).S=0;     
                                agent(SuscepN(i)).IID = InfectN(j);
                                [stateI] = TIV(model,V0(SuscepN(i)));                                
                                agent(SuscepN(i)).TIVs = stateI;                                                    
                                RecoverT = recover_fun(stateI,tspan,params.FoC);
                                agent(SuscepN(i)).RecoverT=tt+RecoverT;
                                stateVL= stateI(1:single(RecoverT/params.IntT)+1,3)';
                                agent(SuscepN(i)).betta(single(tt/params.IntT)+1:single((tt+RecoverT)/params.IntT)+1) = Get_beta(stateVL);
                                agent(SuscepN(i)).taulist = 0:TimeStep:RecoverT;
                                taulist = agent(SuscepN(i)).taulist;
                                agent(SuscepN(i)).Ptest(single(tt/params.IntT)+1:single(tt/params.IntT)+length(taulist)) =...
                                    Prob_get_test(agent(SuscepN(i)).taulist,ModelType,Ph,Pl,TimeStep);
                                agent(SuscepN(i)).Pposit(single(tt/params.IntT)+1:single(tt/params.IntT)+length(taulist)) = ...
                                    Prob_testPosit_VL(agent(SuscepN(i)).TIVs(1:length(taulist),3),minP,range);
                                InfectN=[InfectN SuscepN(i)];
                                SuscepPeople = [SuscepPeople i];                                                             
                            end
                        end
                        SuscepN(SuscepPeople)=[]; % Remove the susceptiable people who have infected
                    end                    
                end
    
            else
                agent(InfectN(j)).I=0;
                agent(InfectN(j)).R=1;
                InfectPeople = [InfectPeople j];
                agent(InfectN(j)).betta(single(tt/params.IntT)+1)=0;
                
            end
        
        end
        InfectN(InfectPeople)=[]; % Remove the infected people who have recovered
        InfectIdtau{single(tt/params.IntT)+1}=InfectN;

    end
    states= [states; sum([agent(:).S]) sum([agent(:).I]) sum([agent(:).R])];
    Ytime=[Ytime tt];
end


end

function [RecoverT] = recover_fun(stateI,tspan,FoC)

    VL_state = stateI(:,end);
    NoT = sum(VL_state>=FoC);
    RecoverT = tspan(NoT);

end

function [isolated] = isolate_policy(OnlyH,HwithL)
    if OnlyH==1&HwithL == 0
        isolated = 1;
    elseif OnlyH==1&HwithL == 1
        isolated = 2;
    elseif OnlyH==0&HwithL == 1
        isolated = 3;
    end

end

