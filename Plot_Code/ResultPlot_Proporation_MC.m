% FLu
%ColorG4 = ["#ffffcc","#edf8b1","#2c7fb8","#253494"];
%ColorG4Pro = ["#fef0d9","#fdd49e","#e34a33","#b30000"];

ColorG4 = ['#dfc27d',"#8c510a","#2c7fb8","#253494"];
ColorG4Pro = ['#80cdc1','#01665e',"#f4a582","#b2182b"];

for i = 1:20
    [a,b] = size(Beta_BeginB);
    for j = 1:a
        Beta_IHP(i,j) = length(Beta_BeginB{j,i});
        Beta_ILP(i,j) = length(Beta_BeginS{j,i});
        Beta_DHP(i,j) = length(Beta_RecoverB{j,i});
        Beta_DLP(i,j) = length(Beta_RecoverS{j,i});
    end
    Total(i,:)= Beta_IHP(i,:)+Beta_ILP(i,:)+Beta_DHP(i,:)+Beta_DLP(i,:);
end

PIH = Beta_IHP./Total;
PIL = Beta_ILP./Total;
PDH = Beta_DHP./Total;
PDL = Beta_DLP./Total;


PIHstd = nanstd(PIH);
PILstd = nanstd(PIL);
PDHstd = nanstd(PDH);
PDLstd = nanstd(PDL);

time = 0:1/24:100;
X_plot = [time, fliplr(time)];
PIL_95 = [mean(PIL,'omitnan')-1.96*PILstd, fliplr(mean(PIL,'omitnan')+1.96*PILstd)];
PDL_95 = [mean(PDL,'omitnan')-1.96*PDLstd, fliplr(mean(PDL,'omitnan')+1.96*PDLstd)];
PIH_95 = [mean(PIH,'omitnan')-1.96*PIHstd, fliplr(mean(PIH,'omitnan')+1.96*PIHstd)];
PDH_95 = [mean(PDH,'omitnan')-1.96*PDHstd, fliplr(mean(PDH,'omitnan')+1.96*PDHstd)];
PIH_95(isnan(PIH_95))=0;
PIL_95(isnan(PIL_95))=0;
PDH_95(isnan(PDH_95))=0;
PDL_95(isnan(PDL_95))=0;


figure(1)
y1 = plot(time,mean(PIL,'omitnan'),'LineWidth',2);
y1.Color = ColorG4(1);
hold on
y1_95 = fill(X_plot,PIL_95,1,'facecolor',ColorG4(1),'edgecolor','none','facealpha',0.3);

y2 = plot(time,mean(PDL,'omitnan'),'LineWidth',2);
y2.Color = ColorG4(2);
y2_95 = fill(X_plot,PDL_95,1,'facecolor',ColorG4(2),'edgecolor','none','facealpha',0.3);

y3 = plot(time,mean(PIH,'omitnan'),'LineWidth',2);
y3.Color = ColorG4(3);
y3_95 = fill(X_plot,PIH_95,1,'facecolor',ColorG4(3),'edgecolor','none','facealpha',0.3);

y4 = plot(time,mean(PDH,'omitnan'),'LineWidth',2);
y4.Color = ColorG4(4);
y4_95 = fill(X_plot,PDH_95,1,'facecolor',ColorG4(4),'edgecolor','none','facealpha',0.3);

xlim([0 60])
ylim([0 1])
xlabel("Time (days)",'fontsize',22)
ylabel("Proportion of infectious population",'fontsize',22)


for j = 1:20
    peak=max(find(states_MC{j}(:,2)==max(states_MC{j}(:,2))));
    peaktime(j) = peak;
end
peakday = ceil(peaktime/24);

lgd1 = legend([y1,y2,y3,y4],'low, increasing','low, decreasing',...
    'high, increasing','high, decreasing','location','northwest','fontsize',12);
lgd1.Title.String = "Viral load";
lgd1.Title.FontSize = 16;
hold on
Q = quantile(peakday,[0.025,0.5,0.975]);
xline(Q,'--','HandleVisibility','off','LineWidth',1.2);
xline(Q(2),'-','HandleVisibility','off','LineWidth',1.2);

box off
legend('boxoff')
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
print(fig,'Pro_MC.pdf','-dpdf')
print(fig,'Pro_MC.eps','-depsc')
print(fig,'Pro_MC.tif','-dtiff')

for i = 1:20
    [a,b] = size(Beta_BeginB);
    for j = 1:a
        Beta_IH(i,j) = sum(Beta_BeginB{j,i});
        Beta_IL(i,j) = sum(Beta_BeginS{j,i});
        Beta_DH(i,j) = sum(Beta_RecoverB{j,i});
        Beta_DL(i,j) = sum(Beta_RecoverS{j,i});
    end
    TotalFoI(i,:)= Beta_IH(i,:)+Beta_IL(i,:)+Beta_DH(i,:)+Beta_DL(i,:);
end

PIHF = Beta_IH./TotalFoI;
PILF = Beta_IL./TotalFoI;
PDHF = Beta_DH./TotalFoI;
PDLF = Beta_DL./TotalFoI;

PIHFstd = nanstd(PIHF);
PILFstd = nanstd(PILF);
PDHFstd = nanstd(PDHF);
PDLFstd = nanstd(PDLF);

PILF_95 = [mean(PILF,'omitnan')-1.96*PILFstd, fliplr(mean(PILF,'omitnan')+1.96*PILFstd)];
PDLF_95 = [mean(PDLF,'omitnan')-1.96*PDLFstd, fliplr(mean(PDLF,'omitnan')+1.96*PDLFstd)];
PIHF_95 = [mean(PIHF,'omitnan')-1.96*PIHFstd, fliplr(mean(PIHF,'omitnan')+1.96*PIHFstd)];
PDHF_95 = [mean(PDHF,'omitnan')-1.96*PDHFstd, fliplr(mean(PDHF,'omitnan')+1.96*PDHFstd)];
PIHF_95(isnan(PIHF_95))=0;
PILF_95(isnan(PILF_95))=0;
PDHF_95(isnan(PDHF_95))=0;
PDLF_95(isnan(PDLF_95))=0;


figure(2)
y1 = plot(time,mean(PILF,'omitnan'),'LineWidth',2);
y1.Color = ColorG4Pro(1);
hold on
y1_95 = fill(X_plot,PILF_95,1,'facecolor',ColorG4Pro(1),'edgecolor','none','facealpha',0.3);

y2 = plot(time,mean(PDLF,'omitnan'),'LineWidth',2);
y2.Color = ColorG4Pro(2);
y2_95 = fill(X_plot,PDLF_95,1,'facecolor',ColorG4Pro(2),'edgecolor','none','facealpha',0.3);

y3 = plot(time,mean(PIHF,'omitnan'),'LineWidth',2);
y3.Color = ColorG4Pro(3);
y3_95 = fill(X_plot,PIHF_95,1,'facecolor',ColorG4Pro(3),'edgecolor','none','facealpha',0.3);

y4 = plot(time,mean(PDHF,'omitnan'),'LineWidth',2);
y4.Color = ColorG4Pro(4);
y4_95 = fill(X_plot,PDHF_95,1,'facecolor',ColorG4Pro(4),'edgecolor','none','facealpha',0.3);

xlim([0 60])
ylim([0 1])
xlabel("Time (days)",'fontsize',22)
ylabel("Proportion of force of infection",'fontsize',22)
lgd2= legend([y1,y2,y3,y4],'low, increasing','low, decreasing',...
    'high, increasing','high, decreasing','location','northwest','fontsize',12);
lgd2.Title.String = "Viral load";
lgd2.Title.FontSize = 16;
fig=gcf;
set(gca,'FontSize',18)
hold on
xline(Q,'--','HandleVisibility','off','LineWidth',1.2);
xline(Q(2),'-','HandleVisibility','off','LineWidth',1.2);


box off
legend('boxoff')
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
print(fig,'ProGrouphl_FluMC.pdf','-dpdf')
print(fig,'ProGrouphl_FluMC.eps','-depsc')
print(fig,'ProGrouphl_FluMC.tif','-dtiff')



PH = Beta_IHP./Total+Beta_DHP./Total;
PL = Beta_ILP./Total+Beta_DLP./Total;

PHstd = nanstd(PH);
PLstd = nanstd(PL);

PL_95 = [mean(PL,'omitnan')-1.96*PLstd, fliplr(mean(PL,'omitnan')+1.96*PLstd)];
PH_95 = [mean(PH,'omitnan')-1.96*PHstd, fliplr(mean(PH,'omitnan')+1.96*PHstd)];
PH_95(isnan(PH_95))=0;
PL_95(isnan(PL_95))=0;

figure(3)
y1 = plot(time,mean(PL,'omitnan'),'LineWidth',2);
y1.Color = ColorG4(2);
hold on
y1_95 = fill(X_plot,PL_95,1,'facecolor',ColorG4(2),'edgecolor','none','facealpha',0.3);

y2 = plot(time,mean(PH,'omitnan'),'LineWidth',2);
y2.Color = ColorG4(4);
y2_95 = fill(X_plot,PH_95,1,'facecolor',ColorG4(4),'edgecolor','none','facealpha',0.3);
xlim([0 60])

ylim([0 1])
xlabel("Time (days)",'fontsize',22)
ylabel("Proportion of infectious population",'fontsize',22)
lgd3= legend([y1,y2],'low','high','location','northwest','fontsize',12);
lgd3.Title.String = "Viral load";
lgd3.Title.FontSize = 16;
hold on
xline(Q,'--','HandleVisibility','off','LineWidth',1.2);
xline(Q(2),'-','HandleVisibility','off','LineWidth',1.2);
box off
legend('boxoff')
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
print(fig,'ProGroup_FluMC.pdf','-dpdf')



PHF = Beta_IH./TotalFoI+Beta_DH./TotalFoI;
PLF = Beta_IL./TotalFoI+Beta_DL./TotalFoI;

PHFstd = nanstd(PHF);
PLFstd = nanstd(PLF);

PLF_95 = [mean(PLF,'omitnan')-1.96*PLFstd, fliplr(mean(PLF,'omitnan')+1.96*PLFstd)];
PHF_95 = [mean(PHF,'omitnan')-1.96*PHFstd, fliplr(mean(PHF,'omitnan')+1.96*PHFstd)];
PHF_95(isnan(PHF_95))=0;
PLF_95(isnan(PLF_95))=0;

figure(4)
y1 = plot(time,mean(PLF,'omitnan'),'LineWidth',2);
y1.Color = ColorG4Pro(2);
hold on
y1_95 = fill(X_plot,PLF_95,1,'facecolor',ColorG4Pro(2),'edgecolor','none','facealpha',0.3);

y2 = plot(time,mean(PHF,'omitnan'),'LineWidth',2);
y2.Color = ColorG4Pro(4);
y2_95 = fill(X_plot,PHF_95,1,'facecolor',ColorG4Pro(4),'edgecolor','none','facealpha',0.3);
xlim([0 60])
ylim([0 1])
xlabel("Time (days)",'fontsize',22)
ylabel("Proportion of force of infection",'fontsize',22)
lgd4= legend([y1,y2],'low','high','location','west','fontsize',12);
lgd4.Title.String = "Viral load";
lgd4.Title.FontSize = 16;
fig=gcf;
set(gca,'FontSize',18)
hold on
xline(Q,'--','HandleVisibility','off','LineWidth',1.2);
xline(Q(2),'-','HandleVisibility','off','LineWidth',1.2);



box off
legend('boxoff')
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf, 'renderer', 'painters');
print(fig,'ProGroupFoI_FluMC.pdf','-dpdf')

