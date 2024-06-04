close all
% clear

% DataMer = load('CORvisco00355.txt');
% DataMervisco0001 = load('NoGravityLBTest/CORvisco0001.txt');
DataMerDry = load('NoGravityLBTest/CORDry174mm.txt');
% DataMerNoGravity = load('NoGravityLBTest/CORNoGravity.txt');
DataWillet = load('NoGravityLBTest/CORWilletandWilletViscous.txt');
DataLS = load('NoGravityLBTest/CORLS.txt');
DataExpMoi = load('ExpMoi.txt');
DataExpDry = load('ExpDry.txt');

% figure
% plot(DataMer(:,1),DataMer(:,3),'--or');
% hold on
% plot(DataMer(:,1),DataMer(:,2),'--ob');
% hold on
% scatter(DataExpDry(:,1),DataExpDry(:,2),'*k');
% hold on
% scatter(DataExpMoi(:,1),DataExpMoi(:,2),'squarek','filled');
% legend('MercuryDPM Dry','MercuryDPM Moist, viscosity=0.03','MercuryDPM Moist, viscosity=0.001','Exp Dry','Exp Moist');
% xlim([0 2]);
% xlabel('Collision velocity [m/s]');ylabel('Restitution Coefficient [-]');

% figure
% plot(DataMerDry(:,1), DataMerDry(:,2), '--ob');
% hold on
% scatter(DataExpDry(:,1),DataExpDry(:,2),'^k','filled');
% xlim([0 2]);
% ylim([0 1]);
% legend('MercuryDPM Dry','Exp Dry');
% xlabel('Collision velocity [m/s]');ylabel('Restitution Coefficient [-]');

figure
% plot(DataMervisco0001(1:6,1),DataMervisco0001(1:6,2),'--om');
% hold on
% plot(DataMervisco0001(7:end,1),DataMervisco0001(7:end,2),'--og');
% hold on
% plot(DataLS(:,1),DataLS(:,2),'--o');
% hold on
% plot(DataWillet(1:6,1),DataWillet(1:6,2),'--o');
% hold on
plot(DataMerDry(:,1),DataMerDry(:,2),'--d','Color',[0 0.4470 0.7410]);
hold on
plot(DataWillet(7:end,1),DataWillet(7:end,2),'--o','Color',[0 0.4470 0.7410]);
hold on
scatter(DataExpDry(:,1),DataExpDry(:,2),'dk','filled');
scatter(DataExpMoi(1:7,1),DataExpMoi(1:7,2),'ok','filled');
% hold on
% scatter(DataExpMoi(8:end,1),DataExpMoi(8:end,2),'squarek','filled');
legend('Simulation Results, dry particle', 'Simulation Results, moist particle', 'Cr$\ddot{u}$ger et al. data, dry particle', 'Cr$\ddot{u}$ger et al. data, moist particle','Interpreter','Latex','fontsize',10);
xlim([0 2]);
xlabel('$U_{\mathrm{impact}}$ [m/s]','Interpreter','Latex','fontsize',18);ylabel('$e_{\mathrm{collision}}$ [-]','Interpreter','Latex','fontsize',18);
set(gca,'FontSize',14);

%%Varying liquid content
DataVaLiquid = load('varyingLiquid/CORVaryVimVliquid.txt');

figure
ax=axes;
ax.ColorOrder=[0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980; 0.6350 0.0780 0.1840; 0.4940 0.1840 0.5560; 0 0 0];
for i = 1:5
    plot(DataVaLiquid(1+5*(i-1):5+5*(i-1),1)*100,DataVaLiquid(1+5*(i-1):5+5*(i-1),2),'--o');
    hold on
end
legend('U_{impact}=0.2 m/s','U_{impact}=0.3 m/s','U_{impact}=0.5 m/s','U_{impact}=1 m/s','U_{impact}=5 m/s');
xlabel('Liquid content by volume [%]');ylabel('Restitution Coefficient [-]');


figure
Vim = [0.2 0.3 0.5 1 5];
for i = 1:5
    COR_Liquid = DataVaLiquid(i:5:i+20,2);
    plot(Vim,COR_Liquid,'--o');
    hold on
end
legend('Dry','Liquid=1%','Liquid=5%','Liquid=10%','Liquid=20%');
xlabel('Impact velocity [m/s]');ylabel('Restitution Coefficient [-]');
xticks(linspace(0,5,6));
