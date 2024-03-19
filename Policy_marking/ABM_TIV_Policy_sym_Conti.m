function [Ytime,states,InfectIdtau,agent,Isolation]=ABM_TIV_Policy_sym_Conti(model,Ytime,States,InfectIdtau,agent,Isolation,IntStartT,DurationT,OnlyH,HwithL,seed)
rng(seed)
beta_tau=model.beta_tau;
Get_beta=beta_tau{1};
condition=model.condition;
params = model.params;
TimeStep = params.IntT;
paramsT=model.paramsT;
tspan =0:paramsT.IntT:paramsT.maxT;
V0 = 10.^(2*rand(1,params.N));
minP = params.minP;
range = params.changeP;



% simulate all the test time
TimeTest = params.TimeTest;

% Record the person is infected at time \tau

InfectN = InfectIdtau{end};
InfectN(ismember(InfectN,Isolation)) = [];
SuscepN=max(InfectIdtau{end})+1:params.N;

states = States;


        
  
for ii = max(Ytime)/TimeStep+1:params.EndT/TimeStep
    tt=ii*TimeStep;

% Need to change this part (first need to get test
    if tt <= DurationT & tt>=IntStartT
        TestInfectN=InfectN(find([agent(InfectN).symtoms]==1));
        for i = 1:length(TestInfectN)
            if agent(TestInfectN(i)).Ptest(ii) == 1
                agent(TestInfectN(i)).TestTime = [agent(TestInfectN(i)).TestTime [tt;0]];
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
                
                timeInfect = tt-agent(InfectN(j)).InfectT;
                timeInfect = max(find(timeInfect>=tspan));
                VL = agent(InfectN(j)).TIVs(timeInfect,3);

                if VL>params.Recover

                    agent(InfectN(j)).betta(single(tt/TimeStep)+1) = Get_beta(VL);

              
                end
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
                                agent(SuscepN(i)).betta(single(tt/params.IntT)+1) = Get_beta(V0(SuscepN(i)));
                                RecoverT = recover_fun(stateI,tspan);
                                agent(SuscepN(i)).RecoverT=tt+RecoverT;
                                agent(SuscepN(i)).taulist = 0:TimeStep:RecoverT;
                                taulist = agent(SuscepN(i)).taulist;
                                TimeRange = 0:ceil(max(taulist))-1;
                                agent(SuscepN(i)).Ptest(single(tt/params.IntT)+1:single(tt/params.IntT)+length(taulist)) = 0;
                                agent(SuscepN(i)).Ptest(single(tt/params.IntT)+TimeTest(SuscepN(i),1:length(TimeRange))+...
                                    1/TimeStep*TimeRange) = 1;
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

function [RecoverT] = recover_fun(stateI,tspan)

    VL_state = stateI(:,end);
    NoT = sum(VL_state>=0.01);
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

