close all
clear

DataMerDry = load('NoGravityLBTest/CORDry174mm.txt');
DataMerWet = load('NoGravityLBTest/WilletViscous200.txt');
DataExpMoi = load('ExpMoi.txt');
DataExpDry = load('ExpDry.txt');

Vim=linspace(0.01,2,200);
Vre=DataMerWet';
CORwet=Vre./Vim;
index_small=find(CORwet<0.088);
CORwet(index_small)=0;

figure
plot(DataMerDry(:,1),DataMerDry(:,2),'--','Color',[0 0.4470 0.7410],'LineWidth',1.5);
hold on;
plot(Vim,CORwet,'-','Color',[0 0.4470 0.7410],'LineWidth',1.5);
scatter(DataExpDry(:,1),DataExpDry(:,2),'dk','filled');
scatter(DataExpMoi(1:7,1),DataExpMoi(1:7,2),'ok','filled');
legend('Simulation Results, dry particle', 'Simulation Results, moist particle', 'Cr$\ddot{u}$ger et al. data, dry particle', 'Cr$\ddot{u}$ger et al. data, moist particle','Interpreter','Latex','fontsize',10);
xlim([0 2]);
xlabel('$U_{\mathrm{impact}}$ [m/s]','Interpreter','Latex','fontsize',18);ylabel('$e_{\mathrm{collision}}$ [-]','Interpreter','Latex','fontsize',18);
set(gca,'FontSize',14);