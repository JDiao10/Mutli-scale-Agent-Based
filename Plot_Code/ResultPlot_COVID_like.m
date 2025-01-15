ColorG4 = ["#e41a1c","#276419","#762a83","#b8e186"];

StepTime = 1/24;
day = single((size(states_MC_B{1},1)-1)*StepTime);
NumNo=100;
N=20000;

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
% 
% 
% row = size(agent_InfectRecoverTime_H,1);
% agent_InfectT_H = zeros(NumNo,row);
% for j = 1:NumNo
%      for i = 1:row
%         Time  = agent_InfectRecoverTime_H{i,j};
%         if length(Time)>0
%             agent_InfectT_H(j,i) = Time(1);
%         else
%             agent_InfectT_H(j,i)= NaN;
%         end
%      end
% end
% 
% 
% row = size(agent_InfectRecoverTime_B,1);
% agent_InfectT_B = zeros(NumNo,row);
% for j = 1:NumNo
%      for i = 1:row
%         Time  = agent_InfectRecoverTime_B{i,j};
%         if length(Time)>0
%             agent_InfectT_B(j,i) = Time(1);
%         else
%             agent_InfectT_B(j,i)= NaN;
%         end
%      end
% end
% 
% 
% row = size(agent_InfectRecoverTime_SF,1);
% agent_InfectT_SF = zeros(NumNo,row);
% for j = 1:NumNo
%      for i = 1:row
%         Time  = agent_InfectRecoverTime_SF{i,j};
%         if length(Time)>0
%             agent_InfectT_SF(j,i) = Time(1);
%         else
%             agent_InfectT_SF(j,i)= NaN;
%         end
%      end
% end

% Single Test
for i = 1:NumNo
    
    for j = 1:200   
        DailyIncidence(i,j) = sum(agent_InfectT(i,:)>=(j-1) & agent_InfectT(i,:)<j);
        % DailyIncidence_H(i,j) = sum(agent_InfectT_H(i,:)>=(j-1) & agent_InfectT_H(i,:)<j);
        % DailyIncidence_B(i,j) = sum(agent_InfectT_B(i,:)>=(j-1) & agent_InfectT_B(i,:)<j);
        % DailyIncidence_SF(i,j) = sum(agent_InfectT_SF(i,:)>=(j-1) & agent_InfectT_SF(i,:)<j);
    end
end

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


figure(1)
time = 1:200;
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
%h = legend("No isolation","Isolate $V_{h}$ only","Isolate all","Isolate $V_{h}$ only if $t>t_{peak}$",...
%    'Interpreter','latex','Fontsize',18,'location','northeast');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
print(fig,'SingleTestCovidDI.pdf','-dpdf')


figure(2)
time = 0:1/24:200;
y1 = plot(time,N-states_MC_No{i}(1:length(time),1),'linewidth',2);
y1.Color = ColorG4(1);

hold on
y2 = plot(time,N-states_MC_H{i}(1:length(time),1),'LineWidth',2);
y2.Color = ColorG4(3);

y4= plot(time,N-states_MC_SF{i}(1:length(time),1),':','LineWidth',2);
y4.Color = ColorG4(4);

y3= plot(time,N-states_MC_B{i}(1:length(time),1),'LineWidth',2);
y3.Color = ColorG4(2);

xlabel("Time (days)",'FontSize',22)
ylabel("Cumulative infections",'FontSize',22)
%h = legend([y1,y2,y3],"No Isolation","Isolate $V_{h}$ only","Isolate all",'Interpreter','latex','Fontsize',18);
%h = legend("No isolation","Isolate $V_{h}$ only","Isolate all","Isolate $V_{h}$ only" + ...
%    "if $t>t_{peak}$",'Interpreter','latex','Fontsize',18);
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

print(fig,'SingleTestCovid.pdf','-dpdf')


s1.Position(1)=0.08;
s1.Position(3) = 0.3;
s1.Position(2) = 0.12;
h.Position(1) = 0.85;
h.Position(3) = 0.07;
s1.Position(4)=0.8;
s2.Position(4)= s1.Position(4);
s2.Position(2)=s1.Position(2);
s2.Position(1)=0.5;
s2.Position(3) = 0.3;

h.Position=[0.816003617568755,0.662449197169862,0.13799276486249,0.227213238106399];

set(gcf, 'PaperSize', [14 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 14 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [14 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 14 6]);
set(gcf, 'renderer', 'painters');

print(fig,'SingleTestCovid.pdf','-dpdf')


%Multiple Tests
%ColorG4=["#e41a1c","#377eb8","#984ea3","#4daf4a"];
NumNo=100;
for i = 1:NumNo
    
    for j = 1:200   
        DailyIncidence(i,j) = sum(agent_InfectT(i,:)>=(j-1) & agent_InfectT(i,:)<j);
        DailyIncidence_H(i,j) = sum(agent_InfectT_H(i,:)>=(j-1) & agent_InfectT_H(i,:)<j);
        DailyIncidence_B(i,j) = sum(agent_InfectT_B(i,:)>=(j-1) & agent_InfectT_B(i,:)<j);
        DailyIncidence_SF(i,j) = sum(agent_InfectT_SF(i,:)>=(j-1) & agent_InfectT_SF(i,:)<j);
    end
end

figure(3)
i = 2;
time = 1:200;
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

h = legend([y1,y3,y2,y4],"No isolation","Isolate all","Isolate high","Adaptive_{0}",...
    'Fontsize',18,'location','northeast');
%h = legend("No isolation","Isolate $V_{h}$ only","Isolate all","Isolate $V_{h}$ only if $t>t_{peak}$",...
%    'Interpreter','latex','Fontsize',18,'location','Northeast');
legend("boxoff")
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
box off
print(fig,'MultiTestCovidDI.pdf','-dpdf')

figure(4)
N = 20000;
time = 0:1/24:200;
y1 = plot(time,N-states_MC_No{i}(1:length(time),1),'linewidth',2);
y1.Color = ColorG4(1);

hold on
y2 = plot(time,N-states_MC_H{i}(1:length(time),1),'LineWidth',2);
y2.Color = ColorG4(3);

y4= plot(time,N-states_MC_SF{i}(1:length(time),1),':','LineWidth',2);
y4.Color = ColorG4(4);

y3= plot(time,N-states_MC_B{i}(1:length(time),1),'LineWidth',2);
y3.Color = ColorG4(2);



xlabel("Time (days)",'FontSize',22)
ylabel("Cumulative infections",'FontSize',22)
%h = legend([y1,y2,y3],"No Isolation","Isolate $V_{h}$ only","Isolate all",'Interpreter','latex','Fontsize',18);
%h = legend("No isolation","Isolate $V_{h}$ only","Isolate all","Isolate $V_{h}$ only" + ...
%    "if $t>t_{peak}$",'Interpreter','latex','Fontsize',18);
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

print(fig,'MultipleTestCovid.pdf','-dpdf')



% effectiveness
for i = 1:100
    e_HFinal(i) = (states_MC_H{i}(end,1)-states_MC_No{i}(end,1))/CummulativeIso_H(i);
    e_BFinal(i) = (states_MC_B{i}(end,1)-states_MC_No{i}(end,1))/CummulativeIso_B(i);
    e_SFFinal(i) = (states_MC_SF{i}(end,1)-states_MC_No{i}(end,1))/CummulativeIso_SF(i);
end

save("EffectivenessMTCovid.mat","e_HFinal","e_BFinal","e_SFFinal")

save("EffectivenessSTCovid.mat","e_HFinal","e_BFinal","e_SFFinal")

figure(4)
N = sum(states(1,:));
time = 0:1/24:200;
y1 = plot(time,N-states(:,1),'linewidth',2);
y1.Color = ColorG4(1);

hold on
y2 = plot(time,N-states_MTH(:,1),'LineWidth',2);
y2.Color = ColorG4(3);


y3= plot(time,N-states_MTB(:,1),'LineWidth',2);
y3.Color = ColorG4(4);

y4= plot(time,N-states_MTSF(:,1),':','LineWidth',2);
y4.Color = ColorG4(4);

xlabel("Time (days)",'FontSize',22)
ylabel("Accumulated cases ($I$)",'Interpreter','latex','FontSize',22)
%h = legend([y1,y2,y3],"No Isolation","Isolate $V_{h}$ only","Isolate all",'Interpreter','latex','Fontsize',18);

fig=gcf;
set(gca,'FontSize',18)

s1.Position(1)=0.08;
s1.Position(3) = 0.3;
s1.Position(2) = 0.12;
h.Position(1) = 0.85;
h.Position(3) = 0.07;
s1.Position(4)=0.8;
s2.Position(4)= s1.Position(4);
s2.Position(2)=s1.Position(2);
s2.Position(1)=0.5;
s2.Position(3) = 0.3;

h.Position=[0.816003617568755,0.662449197169862,0.13799276486249,0.227213238106399];

set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');

print(fig,'MultiTestCovid.pdf','-dpdf')



% asym

% load("MC_Covid_Data.mat")
% load("MC_Covid_Data_PolicyB_MT_sym.mat")
% load("MC_Covid_Data_PolicyH_MT_sym.mat")
% load("MC_Policydata_Covid_StopL_sym_MT.mat")
% 
% ColorG4=["#e41a1c","#377eb8","#984ea3","#4daf4a"];
% 
% StepTime = 1/24;
% day = single((size(states_MC_B{1},1)-1)*StepTime);
% 
% row = size(agent_InfectRecoverTime,1);
% agent_InfectT = zeros(20,row);
% for j = 1:20
%      for i = 1:row
%         Time  = agent_InfectRecoverTime{i,j};
%         if length(Time)>0
%             agent_InfectT(j,i) = Time(1);
%         else
%             agent_InfectT(j,i)= NaN;
%         end
%      end
% end
% 
% 
% row = size(agent_InfectRecoverTime_H,1);
% agent_InfectT_H = zeros(20,row);
% for j = 1:20
%      for i = 1:row
%         Time  = agent_InfectRecoverTime_H{i,j};
%         if length(Time)>0
%             agent_InfectT_H(j,i) = Time(1);
%         else
%             agent_InfectT_H(j,i)= NaN;
%         end
%      end
% end
% 
% 
% row = size(agent_InfectRecoverTime_B,1);
% agent_InfectT_B = zeros(20,row);
% for j = 1:20
%      for i = 1:row
%         Time  = agent_InfectRecoverTime_B{i,j};
%         if length(Time)>0
%             agent_InfectT_B(j,i) = Time(1);
%         else
%             agent_InfectT_B(j,i)= NaN;
%         end
%      end
% end
% 
% 
% row = size(agent_InfectRecoverTime_SF,1);
% agent_InfectT_SF = zeros(20,row);
% for j = 1:20
%      for i = 1:row
%         Time  = agent_InfectRecoverTime_SF{i,j};
%         if length(Time)>0
%             agent_InfectT_SF(j,i) = Time(1);
%         else
%             agent_InfectT_SF(j,i)= NaN;
%         end
%      end
% end
% 
% % Daily test
% for i = 1:20 
% 
%     for j = 1:day    
%         DailyIncidence(i,j) = sum(agent_InfectT(i,:)>=(j-1) & agent_InfectT(i,:)<j);
%         DailyIncidence_H(i,j) = sum(agent_InfectT_H(i,:)>=(j-1) & agent_InfectT_H(i,:)<j);
%         DailyIncidence_B(i,j) = sum(agent_InfectT_B(i,:)>=(j-1) & agent_InfectT_B(i,:)<j);
%         DailyIncidence_SF(i,j) = sum(agent_InfectT_SF(i,:)>=(j-1) & agent_InfectT_SF(i,:)<j);
%     end
% end
% 
% Final_size = [];
% Final_size_H = [];
% Final_size_B = [];
% Final_size_SF = [];
% for i = 1:20
%     state = states_MC{i};
%     state_I(i,:) = state(:,2)';
%     state_H = states_MC_H{i};
%     state_B = states_MC_B{i};
%     state_SF = states_MC_SF{i};
%     Final_size(i,:) = (sum(state,2)-state(:,1))';
%     Final_size_H(i,:) = (sum(state_H,2)-state_H(:,1))';
%     Final_size_B(i,:) = (sum(state_B,2)-state_B(:,1))';
%     Final_size_SF(i,:) = (sum(state_SF,2)-state_SF(:,1))';    
% end

ColorG4=["#e41a1c","#377eb8","#984ea3","#4daf4a"];
NumNo=100;
for i = 1:NumNo
    
    for j = 1:200   
        DailyIncidence(i,j) = sum(agent_InfectT(i,:)>=(j-1) & agent_InfectT(i,:)<j);
        DailyIncidence_H(i,j) = sum(agent_InfectT_H(i,:)>=(j-1) & agent_InfectT_H(i,:)<j);
        DailyIncidence_B(i,j) = sum(agent_InfectT_B(i,:)>=(j-1) & agent_InfectT_B(i,:)<j);
        DailyIncidence_SF(i,j) = sum(agent_InfectT_SF(i,:)>=(j-1) & agent_InfectT_SF(i,:)<j);
    end
end
figure(5)

time = 1:200;
y1 = plot(time,DailyIncidence(i,:),'linewidth',2);
y1.Color = ColorG4(1);

hold on
y2 = plot(time,DailyIncidence_H(i,:),'LineWidth',2);
y2.Color = ColorG4(3);


y3= plot(time,DailyIncidence_B(i,:),'LineWidth',2);
y3.Color = ColorG4(4);

y4= plot(time,DailyIncidence_SF(i,:),':','LineWidth',2);
y4.Color = ColorG4(4);

xlabel("Time (days)",'FontSize',22)
ylabel("Daily incidence",'FontSize',22)
fig=gcf;
set(gca,'FontSize',18)
%h = legend("No isolation","Isolate $V_{h}$ only","Isolate all","Isolate $V_{h}$ only" + ...
%    "if $t>t_{peak}$",'Interpreter','latex','Fontsize',18,'location','northeast');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
xlim([0 150])
print(fig,'DailyTestCovidDI.pdf','-dpdf')



figure(6)
N = sum(states(1,:));
time = 0:1/24:150;
y1 = plot(time,N-states(1:time(end)*24+1,1),'linewidth',2);
y1.Color = ColorG4(1);

hold on
y2 = plot(time,N-states_MTH(1:time(end)*24+1,1),'LineWidth',2);
y2.Color = ColorG4(3);


y3= plot(time,N-states_MTB(1:time(end)*24+1,1),'LineWidth',2);
y3.Color = ColorG4(4);

y4= plot(time,N-states_MTSF(1:time(end)*24+1,1),':','LineWidth',2);
y4.Color = ColorG4(4);

xlabel("Time (days)",'FontSize',22)
ylabel("Accumulated cases ($I$)",'Interpreter','latex','FontSize',22)
%h = legend([y1,y2,y3],"No Isolation","Isolate $V_{h}$ only","Isolate all",'Interpreter','latex','Fontsize',18);

fig=gcf;
set(gca,'FontSize',18)

% s1.Position(1)=0.08;
% s1.Position(3) = 0.3;
% s1.Position(2) = 0.12;
% h.Position(1) = 0.85;
% h.Position(3) = 0.07;
% s1.Position(4)=0.8;
% s2.Position(4)= s1.Position(4);
% s2.Position(2)=s1.Position(2);
% s2.Position(1)=0.5;
% s2.Position(3) = 0.3;
% 
% h.Position=[0.816003617568755,0.662449197169862,0.13799276486249,0.227213238106399];

set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');

print(fig,'DailyTestCovid.pdf','-dpdf')
