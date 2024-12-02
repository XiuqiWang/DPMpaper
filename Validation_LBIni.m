close all
clear
for i=1:6
    eval(['data',num2str(i),'=readMercuryCG("..\WilletViscousTransport\stat\LBIniTP\S00',num2str(i),'DryLBIni.stat");']);
end

D = 0.00025;
dz = 0.00025;%same as h in CG
fac_dim=sqrt((2650/1.225-1)*9.81*D^3)*2650;
Qi = 0.0389777166480000;
Lp=[0.0137916418764536;0.109523873723629;0.408649916341856;0.553683387242477;0.650343656115479;0.762625299646373];
Lpnon = Lp/(2650/1.225*D);
Shields=linspace(0.01,0.06,6);
%Selmani peak position
Lpsel = [0.5 0.8 1.1 1.25 1.5];
Lpsel2 = [0.8 1 1.22 1.5];
Lselnon = [2 3.1 3.5 3.9 4.5];
Lpselnon= Lpsel/(2650/1.225*0.0002);
Lpselnon2= Lpsel2/(2650/1.225*0.0002);
ustarsel = [6 7 8 9 10]*0.13-0.32;
Shieldssel = ustarsel.^2*1.225/(2650-1.225)/9.81/0.0002;

for i=1:6
   eval(['[Q', num2str(i),', Upx', num2str(i),', C', num2str(i),', H', num2str(i),'] = CharTester(data',num2str(i),', ',num2str(i),', D, dz);']); 
   eval(['Qs(',num2str(i),',1)=mean(Q',num2str(i),'(2.5/5*end:end));']); 
   eval(['Qs_std(',num2str(i),',1)=std(Q',num2str(i),'(2.5/5*end:end));']); 
end

%Derive the Q-x data
tx = [0, data1.t(1,1:end)-1];
x=cell(6,5);Qx=cell(6,5);
for i = 1:6
        x{i} = zeros(1, length(Q1)+1); 
        Qx{i} = zeros(1, length(Q1)+1); 
end
for i=1:6
        eval(['Upx_int=[0, Upx',num2str(i),'];']);
        for k=2:length(Q1)+1
            eval(['x{',num2str(i),'}(',num2str(k),')=trapz(tx(1:',num2str(k),'), Upx_int(1:',num2str(k),'));']);
        end
    eval(['Qx{',num2str(i),'}=[Qi, fac_dim*Q',num2str(i),'];']);
end

% %exp data from Creyssels et al.
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
%dry Q-t
t = data1.t(1,:)-1;
figure
subplot(1,2,1)
for i=1:6
    eval(['plot([0, t], [Qi, fac_dim*Q',num2str(7-i),']);']); 
    hold on
end
ylim([0 0.7]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02','$\tilde{\Theta}$=0.01','Interpreter','Latex', 'fontsize',12)
xlabel('$t$ [s]','Interpreter','Latex','fontsize',12);ylabel('$Q$ [kg/m/s]','Interpreter','Latex','fontsize',12);
set(gca,'FontSize',14);
subplot(1,2,2)
errorbar(Shields, Qs, Qs_std, 'or');
hold on;
scatter(Screy(1:end-1), Qcrey(1:end-1),'diamondk')
scatter(SIver(1:end-4), QIver(1:end-4),'*k')
plot(S, Q_B, 'k-.');
plot(S, Q_VR,'k');
plot(x_fit, y_fit, '--r');
text(0.05, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 12);
legend('Simulation results','Creyssels et al. (2009)','Iversen & Rasmussen (1999)','Bagnold (1936)','van Rijn & Strypsteen (2020)','fontsize',10);
xlabel('$\tilde{\Theta}$ [-]','Interpreter','Latex','fontsize',12);
ylabel('$\tilde{Q}_{\mathrm{steady}}$ [-]','Interpreter','Latex','fontsize',12);
box on
grid on
set(gca,'FontSize',14);
ylim([-0.0001 0.04]);

% linearize the x and calculate the corresponding Qx through linear
% interpolation
xlinear=cell(6,1);Qxlinear=cell(6,1);
for i=1:6
        [xnew,id_constantend]=omit_constant_end(x{i});
        Qx_nonconstant=Qx{i}(1:id_constantend);
        eval(['xlinear{',num2str(i),'}=linspace(xnew(1), xnew(end), xnew(end)/0.01);']);
        eval(['Qxlinear{',num2str(i),'}=interp1(xnew, Qx_nonconstant, xlinear{',num2str(i),'});']);
end

%%fitting the Q-x data with the damped oscillator equation in dry condition
initial_guess=[8, pi/500, 0.5; 5, -pi/500, 0.8; 3.5, -pi/100, 0.3; 3.5, -pi/150, 0.3; 3, -pi/100, 0.75;5, -pi/34, 0.8];
lamda=zeros(1,6);phi0=lamda;L=lamda;lamda_error=zeros(6,2);L_error=lamda_error;
for i=1:6
    [lamda(i),phi0(i),L(i),lamda_error(i,:),L_error(i,:)]=Fit(xlinear{i}(1:340),Qxlinear{i}(1:340),Qi,Qs(i),initial_guess(i,:));%[4, -pi/500, 0.8];
end
Q_fitted = zeros(6,500);xideal=linspace(0, 5, 500);
%lamda=initial_guess(:,1);phi0=initial_guess(:,2);L=initial_guess(:,3);
for i=1:6
    Q_fitted(i,:) = Qs(i)+(Qi-Qs(i))*sin(2*pi/lamda(i)*xideal+phi0(i))/sin(phi0(i)) .* exp(-xideal/L(i));
end
Qx_vec = cell2mat(cellfun(@(x) x(1:340), Qxlinear, 'UniformOutput', false));
SS_res = sum((Qx_vec - Q_fitted(:,1:340)).^2); % Sum of squares of residuals
SS_tot = sum((Qx_vec - mean(Q_fitted(:,1:340))).^2); % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

figure
subplot(1,2,1)
box on
for i=1:6
   eval(['plot(x{',num2str(7-i),'}, Qx{',num2str(7-i),'});']); 
   hold on
end
ax = gca;
ax.ColorOrderIndex = 1;
for i=1:6
    plot(xideal, Q_fitted(7-i,:),'--');
end
ylim([0 0.7]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02','$\tilde{\Theta}$=0.01', 'Fitted results', 'Interpreter','Latex', 'fontsize',12)
xlabel('$L$ [m]','Interpreter','Latex','fontsize',12);ylabel('$Q$ [kg/m/s]','Interpreter','Latex','fontsize',12);
set(gca,'FontSize',14);
subplot(1,2,2)
plot(Shields, Lpnon, 'ro', 'LineWidth', 1);
hold on;
plot(Shieldssel, Lpselnon, 'ksquare', 'LineWidth', 1);
plot(Shieldssel(2:end), Lpselnon2, 'kdiamond', 'LineWidth', 1);
text(0.05, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\tilde{\Theta}$ [-]','Interpreter','Latex','fontsize',12);
ylabel('$\tilde{L}_\mathrm{peak}$ [-]','Interpreter','Latex','fontsize',12);
legend('Simulation results, Q_0=0.0390 kg/m/s','Selmani et al., (2018), Q_0=0.0356 kg/m/s','Selmani et al., (2018), Q_0=0.0695 kg/m/s','fontsize',10);
box on
grid on
xlim([0 0.15]);
set(gca,'FontSize',14);

function [a_fitted,b_fitted,c_fitted,a_error,c_error]=Fit(x,Q,Qi,Qs,initial_guess)
options = optimset('TolFun', 1e-8);  % Try lowering or adjusting tolerance
options = optimset('MaxIter', 1000);
% 自定义模型函数
model = @(params, x) Qs+(Qi-Qs)*sin(2*pi/params(1)*x+params(2))/sin(params(2)) .* exp(-x/params(3));
% 使用lsqcurvefit进行拟合
[params_fit, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(model, initial_guess, x, Q);
% Compute confidence intervals for the parameters
ci = nlparci(params_fit, residual, 'jacobian', jacobian);
% 输出拟合结果
a_fitted = params_fit(1);
b_fitted = params_fit(2);
c_fitted = params_fit(3);
a_error = ci(1,:);
c_error = ci(3,:);
end