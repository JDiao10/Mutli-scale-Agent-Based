function [Ytime,states,InfectIdtau,agent]=ABM_TIV(model)

beta_tau=model.beta_tau;
Get_beta=beta_tau{1};
condition=model.condition;
params = model.params;
I0 = params.I0;
TimeStep = params.IntT;
paramsT=model.paramsT;
tspan =0:params.IntT:params.T;
%V0=repmat(9.3*10^-2,params.N,1);
V0 = 10.^(2*rand(1,params.N));

% Record the person is infected at time \tau
InfectIdtau = {};

for i = 1:params.N
    agent(i).I = 0;
    agent(i).S = 1;
    agent(i).R = 0;
    agent(i).InfectT=0;
    agent(i).RecoverT=0;
    agent(i).TIVs = {};
    agent(i).IID = [];
    agent(i).betta=zeros(1,length(tspan));
end

InfectN=1:I0;

for i = 1:length(InfectN)
    agent(i).I=1;
    agent(i).S=0;
    [stateI] = TIV(model,V0(i));
    agent(i).betta(1) = Get_beta(V0(i));
    RecoverT = recover_fun(stateI,tspan);
    agent(i).RecoverT = RecoverT;   
    agent(i).TIVs = stateI;
    agent(i).IID = 0;

end

InfectIdtau{1} = InfectN;

states = [sum([agent(:).S]) sum([agent(:).I]) sum([agent(:).R])];
% Record the Susceptiable individual number
SuscepN=I0+1:params.N;
% Record the time
Ytime=0;
        
  
for tt = TimeStep:TimeStep:params.T
    if ~condition(states(end,:))
        InfectPeople = []; % Record the infected people being recover for each time step
        for j = 1:length(InfectN)
            if tt<=agent(InfectN(j)).RecoverT
                
                timeInfect = tt-agent(InfectN(j)).InfectT;
                timeInfect = max(find(timeInfect>=tspan));
                VL = agent(InfectN(j)).TIVs(timeInfect,3);

                if VL>0.01 
                    %if log10(VL)<0 % only for log10 
                    %    VL = 1;
                    %end
                    agent(InfectN(j)).betta(single(tt/params.IntT)+1) = Get_beta(VL);
              
                end
                if length(SuscepN)>0
                    ContactN = round(poissrnd(params.Arrival*TimeStep));
                    Contact = binornd(ContactN,length(SuscepN)/params.N);
                    %for ii = 1:ContactN
                    %    if rand<length(SuscepN)/params.N
                    %        Contact=Contact+1;
                    %    end
                    %end
                    %Contact = round(ContactN*length(SuscepN)/params.N);
                    if Contact>0
                        SuscepPeople=[];
                        for i = 1:Contact
                            
                            p_ill = agent(InfectN(j)).betta(single(tt/params.IntT)+1)/(params.Arrival);
        
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