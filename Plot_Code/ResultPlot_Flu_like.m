ColorG4=["#e41a1c","#377eb8","#984ea3","#4daf4a"];

ColorG4G = ["#e41a1c","#377eb8","#984ea3","#4daf4a"];



% load no isolation data
StepTime = 1/24;
day = single((size(states_MC_B{1},1)-1)*StepTime);
NumNo=100;
row = size(agent_InfectRecoverTime_No,1);
agent_InfectT = zeros(NumNo,row);
for j = 1:NumNo
     for i = 1:row
        Time  = agent_InfectRecoverTime_No{i,j};
        if length(Time)>0
            agent_InfectT(j,i) = Time(1);
        else
            agent_InfectT(j,i)= NaN;
        end
     end
end

% No isolation
for i = 1:NumNo    
    for j = 1:200   
        DailyIncidence(i,j) = sum(agent_InfectT(i,:)>=(j-1) & agent_InfectT(i,:)<j);
    end
end


% load single test's data
% Different isolation strategies

for i = 1:NumNo
    DailyIncidence_H(i,1)=0;
    DailyIncidence_B(i,1)=0;
    DailyIncidence_SF(i,1)=0;
    for j = 2:200
        DailyIncidence_H(i,j)= -(states_MC_H{i}(j*24+1,1)-states_MC_H{i}((j-1)*24+1,1));
        DailyIncidence_B(i,j)= -(states_MC_B{i}(j*24+1,1)-states_MC_B{i}((j-1)*24+1,1));
        DailyIncidence_SF(i,j)= -(states_MC_SF{i}(j*24+1,1)-states_MC_SF{i}((j-1)*24+1,1));
    end
end

Final_size = [];
Final_size_H = [];
Final_size_B = [];
Final_size_SF = [];
for i = 1:NumNo
    state = states_MC_No{i};
    state_I(i,:) = state(:,2)';
    state_H = states_MC_H{i};
    state_B = states_MC_B{i};
    state_SF = states_MC_SF{i};
    Final_size(i,:) = (sum(state,2)-state(:,1))';
    Final_size_H(i,:) = (sum(state_H,2)-state_H(:,1))';
    Final_size_B(i,:) = (sum(state_B,2)-state_B(:,1))';
    Final_size_SF(i,:) = (sum(state_SF,2)-state_SF(:,1))';    
end

% Sinlge test plot (Fig 3 (c)(d))
figure(1)
ColorG4 = ["#e41a1c","#276419","#762a83","#b8e186"];
i=1;
time = 1:size(DailyIncidence,2);
y1 = plot(time,DailyIncidence(i,:),'linewidth',2);
y1.Color = ColorG4(1);

hold on
y2 = plot(time,DailyIncidence_H(i,:),'LineWidth',2);
y2.Color = ColorG4(3);

y4= plot(time,DailyIncidence_SF(i,:),'LineWidth',2);
y4.Color = ColorG4(4);

y3= plot(time,DailyIncidence_B(i,:),'LineWidth',2);
y3.Color = ColorG4(2);

xlabel("Time (days)",'FontSize',22)
ylabel("Daily incidence",'FontSize',22)
fig=gcf;
set(gca,'FontSize',18)

set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
xlim([0 150])
box off
print(fig,'SingleTestFluDI.pdf','-dpdf')
print(fig,'SingleTestFluDI.eps','-depsc')
print(fig,'SingleTestFluDI.tif','-dtiff')



figure(2)
i=1;
time_2 = 0:StepTime:200;
state = [Final_size(i,:) repmat(Final_size(i,end),1,100*24)];
y1 = plot(time_2,state,'linewidth',2);
y1.Color = ColorG4(1);

hold on
stateH = [Final_size_H(i,:) repmat(Final_size_H(i,end),1,0*24)];
y2 = plot(time_2,stateH,'LineWidth',2);
y2.Color = ColorG4(3);



stateSF = [Final_size_SF(i,:) repmat(Final_size_SF(i,end),1,0*24)];
y4= plot(time_2,stateSF,'LineWidth',2);
y4.Color = ColorG4(4);

stateB = [Final_size_B(i,:) repmat(Final_size_B(i,end),1,0*24)];
y3= plot(time_2,stateB,'LineWidth',2);
y3.Color = ColorG4(2);

xlabel("Time (days)",'FontSize',22)
ylabel("Cumulative infections",'FontSize',22)

fig=gcf;
set(gca,'FontSize',18)


set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
ylim([0 20000])
xlim([0 150])
box off
print(fig,'SingleTestFlu.pdf','-dpdf')
print(fig,'SingleTestFlu.eps','-depsc')
print(fig,'SingleTestFlu.tif','-dtiff')

for i = 1:100
    ePerTime_H(i) = (Final_size(i,end)-Final_size_H(i,end))/sum(IsolationTime_H{i});
    ePerTime_B(i) = (Final_size(i,end)-Final_size_B(i,end))/sum(IsolationTime_B{i});
    ePerTime_SF(i) = (Final_size(i,end)-Final_size_SF(i,end))/sum(IsolationTime_SF{i});
end
% Efficiency for Single tests
save("EfficiencyPerDayST.mat",'ePerTime_H','ePerTime_B','ePerTime_SF')

for i = 1:100
    phi_Isolation(i) = 1-(IsolationPer_SF(i,end)-IsolationPer_H(i,end))./(IsolationPer_B(i,end)-IsolationPer_H(i,end));

end
save("Isolation_PerST.mat","IsolationPer_H","IsolationPer_B","IsolationPer_SF")


% Multiple tests
% load multiple tests' data
for i = 1:NumNo
    DailyIncidence_H(i,1)=0;
    DailyIncidence_B(i,1)=0;
    DailyIncidence_SF(i,1)=0;
    for j = 2:200
        DailyIncidence_H(i,j)= -(states_MC_H{i}(j*24+1,1)-states_MC_H{i}((j-1)*24+1,1));
        DailyIncidence_B(i,j)= -(states_MC_B{i}(j*24+1,1)-states_MC_B{i}((j-1)*24+1,1));
        DailyIncidence_SF(i,j)= -(states_MC_SF{i}(j*24+1,1)-states_MC_SF{i}((j-1)*24+1,1));
    end
end

Final_size = [];
Final_size_H = [];
Final_size_B = [];
Final_size_SF = [];
for i = 1:NumNo
    state = states_MC_No{i};
    state_I(i,:) = state(:,2)';
    state_H = states_MC_H{i};
    state_B = states_MC_B{i};
    state_SF = states_MC_SF{i};
    Final_size(i,:) = (sum(state,2)-state(:,1))';
    Final_size_H(i,:) = (sum(state_H,2)-state_H(:,1))';
    Final_size_B(i,:) = (sum(state_B,2)-state_B(:,1))';
    Final_size_SF(i,:) = (sum(state_SF,2)-state_SF(:,1))';    
end
ColorG4=["#e41a1c","#377eb8","#984ea3","#4daf4a"];
ColorG4G = ["#e41a1c","#377eb8","#984ea3","#4daf4a"];

NumNo = 100;

% Mutliple tests plot (Fig 3(a)(b))
figure(3)
i = 12;
time = 1:200;
ColorG4 = ["#e41a1c","#276419","#762a83","#b8e186"];
%ColorG4 = ["#e41a1c","#7f3b08","#984ea3","#fee0b6"];

y1 = plot(time,DailyIncidence(i,:),'linewidth',2);
y1.Color = ColorG4(1);

hold on
y2 = plot(time,DailyIncidence_H(i,:),'LineWidth',2);
y2.Color = ColorG4(3);

y4= plot(time,DailyIncidence_SF(i,:),'LineWidth',2);
y4.Color = ColorG4(4);

y3= plot(time,DailyIncidence_B(i,:),'LineWidth',2);
y3.Color = ColorG4(2);



xlabel("Time (days)",'FontSize',22)
ylabel("Daily incidence",'FontSize',22)

h = legend([y1,y3,y2,y4],"No isolation","Isolate all","Isolate high","Adaptive_{0}",...
    'Fontsize',18,'location','northeast');
legend('boxoff')
fig=gcf;
set(gca,'FontSize',18)
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
xlim([0 150])
box off
print(fig,'DailyIncidenceFluMT_v1.pdf','-dpdf')
print(fig,'DailyIncidenceFluMT_v1.eps','-depsc')
print(fig,'DailyIncidenceFluMT_v1.tif','-dtiff')


figure(4)
i=12;
time_2 = 0:StepTime:250;
state = [Final_size(i,:) repmat(Final_size(i,end),1,150*24)];
y1 = plot(time_2,state,'linewidth',2);
y1.Color = ColorG4(1);

hold on
stateH = [Final_size_H(i,:) repmat(Final_size_H(i,end),1,0*24)];
y2 = plot(time_2,stateH,'LineWidth',2);
y2.Color = ColorG4(3);



stateSF = [Final_size_SF(i,:) repmat(Final_size_SF(i,end),1,0*24)];
y4= plot(time_2,stateSF,'LineWidth',2);
y4.Color = ColorG4(4);

stateB = [Final_size_B(i,:) repmat(Final_size_B(i,end),1,0*24)];
y3= plot(time_2,stateB,'LineWidth',2);
y3.Color = ColorG4(2);

xlabel("Time (days)",'FontSize',22)
ylabel("Cumulative infections",'FontSize',22)

fig=gcf;
set(gca,'FontSize',18)

set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
ylim([0 20000])
xlim([0 150])
box off
print(fig,'MultipleTestFlu.pdf','-dpdf')
print(fig,'MultipleTestFlu.eps','-depsc')
print(fig,'MultipleTestFlu.tif','-dtiff')



for i = 1:100
    ePerTime_H(i) = (Final_size(i,end)-Final_size_H(i,end))/sum(IsolationTime_H{i});
    ePerTime_B(i) = (Final_size(i,end)-Final_size_B(i,end))/sum(IsolationTime_B{i});
    ePerTime_SF(i) = (Final_size(i,end)-Final_size_SF(i,end))/sum(IsolationTime_SF{i});
end
% Efficiency for multiple tests
save("EfficiencyPerDayMT.mat",'ePerTime_H','ePerTime_B','ePerTime_SF')




% Calculate the relative cost


for i = 1:100
    phi_Isolation(i) = 1-(IsolationPer_SF(i,end)-IsolationPer_H(i,end))./(IsolationPer_B(i,end)-IsolationPer_H(i,end));

end
save("Isolation_PerMT.mat","IsolationPer_H","IsolationPer_B","IsolationPer_SF")







% Save 95% quantile of the effectiveness
load('EffectivenessMT.mat')
load('EffectivenessST.mat')

H_95 = quantile(e_HFinal,[0.025 0.975]);
B_95 = quantile(e_BFinal,[0.025 0.975]);
SF_95 = quantile(e_SFFinal,[0.025 0.975]);
e_HFinal95 = e_HFinal(find(e_HFinal>=H_95(1)& e_HFinal<H_95(2)))';
e_BFinal95 = e_BFinal(find(e_BFinal>=B_95(1)& e_BFinal<B_95(2)))';
e_SFFinal95 = e_SFFinal(find(e_SFFinal>=SF_95(1)& e_SFFinal<SF_95(2)))';

save("EffectivenessMT95.mat","e_HFinal95","e_BFinal95","e_SFFinal95")
save("EffectivenessST95.mat","e_HFinal95","e_BFinal95","e_SFFinal95")


% Asymptomatic case

load("Policydata_NoP.mat")
load("Policydata_MTH_sym.mat")
load("Policydata_AllMT_sym.mat")
load("Policydata_MTStopL_sym.mat")

ColorG4=["#e41a1c","#377eb8","#984ea3","#4daf4a"];

for i = 1:NumNo
    
    for j = 1:200   
        DailyIncidence(i,j) = sum(agent_InfectT(i,:)>=(j-1) & agent_InfectT(i,:)<j);
        DailyIncidence_H(i,j) = sum(agent_InfectT_H(i,:)>=(j-1) & agent_InfectT_H(i,:)<j);
        DailyIncidence_B(i,j) = sum(agent_InfectT_B(i,:)>=(j-1) & agent_InfectT_B(i,:)<j);
        DailyIncidence_SF(i,j) = sum(agent_InfectT_SF(i,:)>=(j-1) & agent_InfectT_SF(i,:)<j);
    end
end

figure(5)

i = 1;
time = 1:200;
y1 = plot(time,DailyIncidence(i,:),'linewidth',2);
y1.Color = ColorG4(1);

hold on
y3= plot(time,DailyIncidence_B(i,:),'LineWidth',2);
y3.Color = ColorG4(4);

y2 = plot(time,DailyIncidence_H(i,:),'LineWidth',2);
y2.Color = ColorG4(3);


y4= plot(time,DailyIncidence_SF(i,:),':','LineWidth',2);
y4.Color = ColorG4(4);

xlabel("Time (days)",'FontSize',22)
ylabel("Daily incidence",'FontSize',22)

h = legend("No isolation","Isolate all","Isolate high","$\textrm{Adaptive}_0$",...
    'Interpreter','latex','Fontsize',18,'location','northeast');

xlabel("Time (days)",'FontSize',22)
ylabel("Daily incidence",'FontSize',22)
fig=gcf;
set(gca,'FontSize',18)
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
xlim([0 150])
print(fig,'DailyTesttFluDI.pdf','-dpdf')


Final_size = [];
Final_size_H = [];
Final_size_B = [];
Final_size_SF = [];
for i = 1:NumNo
    state = states_MC_No{i};
    state_I(i,:) = state(:,2)';
    state_H = states_MC_H{i};
    state_B = states_MC_B{i};
    state_SF = states_MC_SF{i};
    Final_size(i,:) = (sum(state,2)-state(:,1))';
    Final_size_H(i,:) = (sum(state_H,2)-state_H(:,1))';
    Final_size_B(i,:) = (sum(state_B,2)-state_B(:,1))';
    Final_size_SF(i,:) = (sum(state_SF,2)-state_SF(:,1))';    
end
figure(6)
i=1;
time_2 = 0:StepTime:200;
state = [Final_size(i,:) repmat(Final_size(i,end),1,100*24)];
y1 = plot(time_2,state,'linewidth',2);
y1.Color = ColorG4(1);

hold on
stateH = [Final_size_H(i,:) repmat(Final_size_H(i,end),1,0*24)];
y2 = plot(time_2,stateH,'LineWidth',2);
y2.Color = ColorG4(3);

stateB = [Final_size_B(i,:) repmat(Final_size_B(i,end),1,0*24)];
y3= plot(time_2,stateB,'LineWidth',2);
y3.Color = ColorG4(4);

stateSF = [Final_size_SF(i,:) repmat(Final_size_SF(i,end),1,0*24)];
y4= plot(time_2,stateSF,':','LineWidth',2);
y4.Color = ColorG4(4);


xlabel("Time (days)",'FontSize',22)
ylabel("Cumulative cases ($I$)",'Interpreter','latex','FontSize',22)

fig=gcf;
set(gca,'FontSize',18)
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
xlim([0 150])
ylim([0 20000])
print(fig,'DailyTesttFlu.pdf','-dpdf')