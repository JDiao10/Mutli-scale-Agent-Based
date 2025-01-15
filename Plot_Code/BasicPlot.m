dt = 1/72;
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
model.paramsT = paramsT;

%ColorG4 = ["#ffffcc","#edf8b1","#2c7fb8","#253494"];
ColorG4 = ['#dfc27d',"#8c510a","#2c7fb8","#253494"];
V0 = 1;
[TIVstate] = TIV(model,V0);
TIVstate(:,3) = TIVstate(:,3)+6000;
TIVtime= 0:360-1;
%pp = plot(TIVstate(1:length(TIVtime),3),'linewidth',2);

yV = TIVstate(1:length(TIVtime),3);
ix = 99;
iy = yV(ix);


ixmax = 126;
iymax = yV(ixmax);


ixmaxb = find(diff(sign(yV-iymax)));
ixmaxb = max(ixmaxb);

xpeak = find(TIVstate(1:length(TIVtime),3)==max(TIVstate(:,3)));
ypeak = yV(xpeak);


lw = 9;

figure(1)
plot(1:TIVtime(end),TIVstate(1:TIVtime(end),3)','Color','k','LineWidth',lw+1)
hold on 
p1 = plot(1:ixmax,TIVstate(1:ixmax,3),'linewidth',lw);
p1.Color = ColorG4(1);

peak = find(TIVstate(:,3)==max(TIVstate(:,3)));
p2 = plot(ixmax:peak,TIVstate(ixmax:peak,3),'LineWidth',lw);
p2.Color = ColorG4(3);


p3 = plot(peak:ixmaxb,TIVstate(peak:ixmaxb,3),'LineWidth',lw);
p3.Color = ColorG4(4);


p4 = plot(ixmaxb+1:TIVtime(end),TIVstate(ixmaxb+1:TIVtime(end),3),'LineWidth',lw);
p4.Color = ColorG4(2);


% annotation('textbox',[0.28,0.1,0.1,0.2],'String','$V_{\ell}^{inc}$',...
%     'Interpreter','latex','FontSize',22,'EdgeColor','none')
% 
% annotation('textbox',[0.33,0.6,0.1,0.2],'String','$V_{h}^{inc}$',...
%     'Interpreter','latex','FontSize',22,'EdgeColor','none')
% 
% annotation('textbox',[0.48,0.6,0.1,0.2],'String','$V_{h}^{dec}$',...
%     'Interpreter','latex','FontSize',22,'EdgeColor','none')
% 
% annotation('textbox',[0.60,0.1,0.1,0.2],'String','$V_{\ell}^{dec}$',...
%     'Interpreter','latex','FontSize',22,'EdgeColor','none')


annotation('textbox',[0.24,0.15,0.1,0.2],'String','low, increasing',...
    'Interpreter','latex','FontSize',22,'EdgeColor','none')

annotation('textbox',[0.28,0.6,0.1,0.2],'String','high, increasing',...
    'Interpreter','latex','FontSize',22,'EdgeColor','none')

annotation('textbox',[0.48,0.6,0.1,0.2],'String','high, decreasing',...
    'Interpreter','latex','FontSize',22,'EdgeColor','none')

annotation('textbox',[0.56,0.15,0.1,0.2],'String','low, decreasing',...
    'Interpreter','latex','FontSize',22,'EdgeColor','none')


plot([0 ixmax], [1 1]*iymax,'--k','linewidth',2)
plot([1 1]*ixmax,[0 1]*iymax,'--k','linewidth',2)
plot(ixmax,iymax,'ok')
text(-1,iymax,'$V^*$','Interpreter','latex','Horiz','right', 'Vert','middle','fontsize',22)
text(ixmax,-2000,'$\tau_{inc}^{*}$','Interpreter','latex', 'Horiz','center', 'Vert','top','fontsize',22)
plot([0 ixmaxb+1], [1 1]*iymax,'--k','linewidth',2)
plot([1 1]*ixmaxb+1,[0 1]*iymax,'--k','linewidth',2)
plot(ixmaxb+1,iymax,'ok')
text(ixmaxb+1,-2000,'$\tau_{dec}^{*}$','Interpreter','latex', 'Horiz','center', 'Vert','top','fontsize',22)
set(gca,'xtick',[],'ytick',[])
labelx = xlabel('\rm{Time} ($\tau$)','Interpreter','latex','fontSize',22);
labely = ylabel('Viral load','fontSize',22);
labely.Position = [-14,iymax];
labelx.Position = [ixmaxb, -30000];

plot([1 1]*xpeak, [0 1]*ypeak,'--k','LineWidth',2)
plot(xpeak,ypeak,'ok')
text(xpeak,-2000,'$\tau_{peak}$','Interpreter','latex', 'Horiz','center', 'Vert','top','fontsize',22)



fig=gcf;
set(gca,'FontSize',18)
set(gcf, 'PaperSize', [16 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 16 8]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [16 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 16 8]);
set(gcf, 'renderer', 'painters');
box off
print(fig,'Basic_v2.pdf','-dpdf')


time = 0:9;
Pro = [repmat(0.2,1,2) repmat(0.8,1,8)];
figure(2)
scatter(time,Pro,'ob','LineWidth',2)
xlabel('\textrm{Infectious time} ($\tau$)','Interpreter','latex','FontSize',18)
ylabel('$P(\textrm{Tested})$','FontSize',18,'Interpreter','latex')
fig=gcf;
set(gca,'FontSize',18)
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 8]);
set(gcf, 'renderer', 'painters');
ylim([0 1])
xlim([0 9])
box off
print(fig,'ProTest_Flu.pdf','-dpdf')


time = 0:21;
Pro = [repmat(0.1,1,3) repmat(0.5,1,19)];
figure(2)
scatter(time,Pro,'ob','LineWidth',2)
xlabel('\textrm{Infectious time} ($\tau$)','Interpreter','latex','FontSize',18)
ylabel('$P(\textrm{Tested})$','FontSize',18,'Interpreter','latex')
fig=gcf;
set(gca,'FontSize',18)
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 8]);
set(gcf, 'renderer', 'painters');
ylim([0 1])
xlim([0 21])
box off
print(fig,'ProTest_Covid.pdf','-dpdf')