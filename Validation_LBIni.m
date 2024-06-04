close all
% clear
% for i=1:6
%     eval(['data',num2str(i),'=readMercuryCG("..\WilletViscousTransport\stat\LBIniTP\S00',num2str(i),'DryLBIni.stat");']);
% end
D = 0.00025;
dz = 0.00025;%same as h in CG
fac_dim=sqrt((2650/1.225-1)*9.81*D^3)*2650;
for i=1:6
   eval(['[Q', num2str(i),', V', num2str(i),', C', num2str(i),', H', num2str(i),'] = CharTester(data',num2str(i),', D, dz);']); 
   eval(['Qs(',num2str(i),',1)=mean(Q',num2str(i),'(2.5/5*end:end));']); 
   eval(['Qs_std(',num2str(i),',1)=std(Q',num2str(i),'(2.5/5*end:end));']); 
end

%dry Q-t
figure
for i=1:6
    eval(['plot(data',num2str(7-i),'.t(1,:)-1, fac_dim*Q',num2str(7-i),');']); 
    hold on
end
ylim([0 0.7]);
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02','$\tilde{\Theta}$=0.01', 'Interpreter','Latex', 'fontsize',12)
xlabel('$t$ [s]','Interpreter','Latex','fontsize',12);ylabel('$Q$ [kg/m/s]','Interpreter','Latex','fontsize',12);


%exp data from Creyssels et al.
Screy = [0.01 0.0225 0.035 0.05 0.0675 0.0975];
Qcrey = [0.005 0.01 0.016 0.024 0.0375 0.06];%dimensionless
%exp data from Iversen&Rasmussen
SIver = [0.01 0.015 0.035 0.03 0.035 0.045 0.045 0.055 0.055 0.0575 0.08 0.0875 0.095 0.098 0.078];
QIver = [0.0001 0.0025 0.0075 0.0125 0.018 0.0125 0.022 0.013 0.0185 0.0375 0.0225 0.028 0.025 0.0275 0.058];

%Empirical formulations
S = linspace(0,0.08,9);
%Bagnold
u_star = sqrt(S*(2650-1.225)*9.81*D/1.225);
u_st = 0.085*sqrt(9.81*D*(2650-1.225)/1.225);

Q_B = 1.8*1.225/9.81*u_star.^3./sqrt((2650/1.225-1)*9.81*D^3)/2650;
%Zingg
Q_Z = 0.83*1.225/9.81*u_star.^3;
%Kawamura
Q_K = 2.78*1.225/9.81*(u_star-u_st).*(u_star+0.234).^2./sqrt((2650/1.225-1)*9.81*D^3)/2650;
%Lettau and Lettau
Q_LL = 6.7*1.225/9.81*u_star.^2.*(u_star-u_st)./sqrt((2650/1.225-1)*9.81*D^3)/2650;
%Hsu
C_H = exp(-0.47+4.97*D)/10000;
Q_H = C_H*(u_star./sqrt(9.81*D)).^3./sqrt((2650/1.225-1)*9.81*D^3)/2650;
%van rijn (1)
Q_VR = 1.5*1.225/9.81*((2650-1.225)*9.81*D/1.225)^1.5*(S.^1.5-0.11^3);

Shields=linspace(0.01,0.06,6);
p = polyfit(Shields, Qs, 1); 
x_fit = linspace(min(Shields), max(Shields), 100); 
y_fit = polyval(p, x_fit); 

%%validation figure for dry limit
figure
errorbar(Shields, Qs, Qs_std, 'or');
hold on
scatter(Screy(1:end-1), Qcrey(1:end-1),'diamondk')
hold on
scatter(SIver(1:end-4), QIver(1:end-4),'*k')
hold on
plot(S, Q_B, 'k-.');
hold on
plot(S, Q_VR,'k');
plot(x_fit, y_fit, '--r');
legend('Simulation results','Creyssels et al. data','Iversen & Rasmussen data','Bagnold model','van Rijn & Strypsteen model','fontsize',10);
xlabel('$\tilde{\Theta}$ [-]','Interpreter','Latex','fontsize',12);
ylabel('$\tilde{Q}_{\mathrm{steady}}$ [-]','Interpreter','Latex','fontsize',12);
box on
grid on
set(gca,'FontSize',14);
ylim([-0.0001 0.04]);

function smoothed_data=smoothdata(Q)
data = Q;
% Define window size for moving average
window_size = 30;
% Apply moving average filter
smoothed_data = movmean(data, window_size);
end