close all
clear
data10=readMercuryCG('stat\LBIniTP\S001DryLBIni.stat');
data11=readMercuryCG('stat\LBIniTP\S001M1LBIni.stat');
data12=readMercuryCG('stat\LBIniTP\S001M5LBIni.stat');
data13=readMercuryCG('stat\LBIniTP\S001M10LBIni.stat');
data14=readMercuryCG('stat\LBIniTP\S001M20LBIni.stat');

data20=readMercuryCG('stat\LBIniTP\S002DryLBIni.stat');
data21=readMercuryCG('stat\LBIniTP\S002M1LBIni.stat');
data22=readMercuryCG('stat\LBIniTP\S002M5LBIni.stat');
data23=readMercuryCG('stat\LBIniTP\S002M10LBIni.stat');
data24=readMercuryCG('stat\LBIniTP\S002M20LBIni.stat');

data30=readMercuryCG('stat\LBIniTP\S003DryLBIni.stat');
data31=readMercuryCG('stat\LBIniTP\S003M1LBIni.stat');
data32=readMercuryCG('stat\LBIniTP\S003M5LBIni.stat');
data33=readMercuryCG('stat\LBIniTP\S003M10LBIni.stat');
data34=readMercuryCG('stat\LBIniTP\S003M20LBIni.stat');
 
data40=readMercuryCG('stat\LBIniTP\S004DryLBIni.stat');
data41=readMercuryCG('stat\LBIniTP\S004M1LBR1Ini.stat');
data42=readMercuryCG('stat\LBIniTP\S004M5LBIni.stat');
data43=readMercuryCG('stat\LBIniTP\S004M10LBIni.stat');
data44=readMercuryCG('stat\LBIniTP\S004M20LBIni.stat');

data50=readMercuryCG('stat\LBIniTP\S005DryLBIni.stat');
data51=readMercuryCG('stat\LBIniTP\S005M1LBIni.stat');
data52=readMercuryCG('stat\LBIniTP\S005M5LBR1Ini.stat');
data53=readMercuryCG('stat\LBIniTP\S005M10LBR1Ini.stat');
data54=readMercuryCG('stat\LBIniTP\S005M20LBIni.stat');

data60=readMercuryCG('stat\LBIniTP\S006DryLBIni.stat');
data61=readMercuryCG('stat\LBIniTP\S006M1LBIni.stat');
data62=readMercuryCG('stat\LBIniTP\S006M5LBIni.stat');
data63=readMercuryCG('stat\LBIniTP\S006M10LBIni.stat');
data64=readMercuryCG('stat\LBIniTP\S006M20LBIni.stat');

D = 0.00025;
dz = 0.00025;%same as h in CG
Liquid=[0 1 5 10 20];
Shields=linspace(0.01,0.06,6);
ustar=sqrt(Shields*(2650-1.225)*9.81*D/1.225);
ustart=sqrt(0.005*(2650-1.225)*9.81*D/1.225);
Qi=4.87221458100000e-07;
dt=0.01;
%different initial particle flux
Qithree=[6.28867960993324e-08,2.96633000380294e-07,4.87221458100000e-07];
Lpeakthree=[0.829199326490001,8.20266023248140,23.9509624289844,40.3470462563418,38.6608653976041;0.709328162292667,4.14300014457050,9.22045260007650,13.1565687250352,14.8297041997901;0.776109244734347,2.02348632879845,7.00235163086276,9.92093851696785,13.9562134004145];
colors = [0 0 0;0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];

for j=1:6
for i=0:4
    eval(['[Q',num2str(j),num2str(i),', Upx',num2str(j), num2str(i), ', Cvx',num2str(j), num2str(i),', H',num2str(j), num2str(i),'] = CharTester(data',num2str(j), num2str(i),',',num2str(j),',D, dz);']); 
end
end

%Qs,Qstd,Cs,Ups
fac_dim=sqrt((2650/1.225-1)*9.81*D^3)*2650;
for i=1:6
     eval(['Qs(',num2str(i),')=fac_dim*mean(Q',num2str(i),'0(2.5/5*end:end));']);
end

for i=1:6
for j=0:4
if i<4
   eval(['Qsm(',num2str(i),',',num2str(j+1),')=fac_dim*mean(Q',num2str(i),num2str(j),'(2.5/5*end:end));']);
   eval(['Cs(',num2str(i),',',num2str(j+1),')=mean(Cvx',num2str(i),num2str(j),'(2.5/5*end:end));']);
   eval(['[Ups(',num2str(i),',',num2str(j+1),'),Upsstd(',num2str(i),',',num2str(j+1),')]=getMeanSTDOfNonZero(Upx',num2str(i),num2str(j),'(2.5/5*end:end));']);
   eval(['[Hs(',num2str(i),',',num2str(j+1),'),Hsstd(',num2str(i),',',num2str(j+1),')]=getMeanSTDOfNonZero(H',num2str(i),num2str(j),'(2.5/5*end:end));']);
   eval(['Qsstd(',num2str(i),',',num2str(j+1),')=fac_dim*std(Q',num2str(i),num2str(j),'(2.5/5*end:end));']);
   eval(['Csstd(',num2str(i),',',num2str(j+1),')=std(Cvx',num2str(i),num2str(j),'(2.5/5*end:end));']);
   else
       eval(['Qsm(',num2str(i),',',num2str(j+1),')=fac_dim*mean(Q',num2str(i),num2str(j),'(3.5/5*end:end));']); 
       eval(['Cs(',num2str(i),',',num2str(j+1),')=mean(Cvx',num2str(i),num2str(j),'(3.5/5*end:end));']);
       eval(['[Ups(',num2str(i),',',num2str(j+1),'),Upsstd(',num2str(i),',',num2str(j+1),')]=getMeanSTDOfNonZero(Upx',num2str(i),num2str(j),'(3.5/5*end:end));']);
       eval(['[Hs(',num2str(i),',',num2str(j+1),'), Hsstd(',num2str(i),',',num2str(j+1),')]=getMeanSTDOfNonZero(H',num2str(i),num2str(j),'(3.5/5*end:end));']);
       eval(['Qsstd(',num2str(i),',',num2str(j+1),')=fac_dim*std(Q',num2str(i),num2str(j),'(3.5/5*end:end));']);
       eval(['Csstd(',num2str(i),',',num2str(j+1),')=std(Cvx',num2str(i),num2str(j),'(3.5/5*end:end));']);
end
end
end

figure
subplot(2,2,1)
for i=1:size(Qsm,1)
errorbar(Liquid,Qsm(length(Shields)-i+1,:),Qsstd(length(Shields)-i+1,:),'o','LineWidth', 1);
hold on
end
plot(Liquid,mean(Qsm),'--','Color', [0.5, 0.5, 0.5]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
ylim([-0.01 0.07]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$Q_{\mathrm{steady}}$ [kg/m/s]','Interpreter','Latex','fontsize',12);
subplot(2,2,2)
for i=1:size(Qsm,1)
errorbar(Liquid,Cs(length(Shields)-i+1,:),Csstd(length(Shields)-i+1,:),'o','LineWidth', 1);
hold on
end
plot(Liquid,mean(Cs),'--','Color', [0.5, 0.5, 0.5]);
text(0.05, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$C_{\mathrm{sal,steady}}$ [kg/m$^2$]','Interpreter','Latex','fontsize',12);
subplot(2,2,3)
for i=1:size(Qsm,1)
errorbar(Liquid,Ups(length(Shields)-i+1,:),Upsstd(length(Shields)-i+1,:),'o','LineWidth', 1);
hold on
end
plot(Liquid,mean(Ups),'--','Color', [0.5, 0.5, 0.5]);
text(0.05, 0.95, '(c)', 'Units', 'normalized', 'FontSize', 12);
ylim([0 3.75]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$U_{\mathrm{sal,steady}}$ [m/s]','Interpreter','Latex','fontsize',12);
subplot(2,2,4)
for i=1:size(Qsm,1)
errorbar(Liquid,Hs(length(Shields)-i+1,:),Hsstd(length(Shields)-i+1,:),'o','LineWidth', 1);
hold on
end
plot(Liquid,mean(Hs),'--','Color', [0.5, 0.5, 0.5]);
text(0.05, 0.95, '(d)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$H_{\mathrm{50,steady}}$ [m]','Interpreter','Latex','fontsize',12);
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02', '$\tilde{\Theta}$=0.01', 'Interpreter','Latex', 'fontsize',10);

%peak and saturated time scales from Q-t
t = data60.t(1,:)-1;
for i=1:6
    for j=0:4
        eval(['[Qtpeak,t_peak, index_peak, t_peak_error]=getIndexPeak(t, fac_dim*Q',num2str(i),num2str(j),');']);
        t_peak_mean(i,j+1)=t_peak;
        t_peak_std(i,j+1)=t_peak_error;
        eval(['Qpeak(',num2str(i),',',num2str(j+1),')=Qtpeak;']);
        %search from index_peak
        eval(['[t_sat,t_sat_error]=getSaturationIndex(t(index_peak:end), fac_dim*Q',num2str(i),num2str(j),'(index_peak:end),Qsm(',num2str(i),'));']);
        t_dep_mean(i,j+1)=t_sat-t_peak;
%         t_sat_std(i,j+1)=t_sat_error;
    end
end

figure
for i=1:6
subplot(3,2,i);
box on
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for j=0:4
   eval(['plot([0, t], [Qi, fac_dim*Q',num2str(i),num2str(j),']);']); 
   hold on
end
ylim([0, 0.7]);xlim([0, 5]);
xlabel('$t$ [s]','Interpreter','Latex','fontsize',12);ylabel('$Q$ [kg/m/s]','Interpreter','Latex','fontsize',12);
eval(['title("$\tilde{\Theta}$=0.0',num2str(i),'","Interpreter","Latex");']);
if i==1
   legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);
end
end

figure
subplot(1,2,1)
for i = 1:6
plot(Liquid, t_peak_mean(7-i,:), 'o', 'LineWidth', 1);
hold on
end
ylim([0 3.5]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$T_{\mathrm{peak}}$ [s]','Interpreter','Latex','fontsize',12);
subplot(1,2,2)
for i = 1:6
plot(Liquid, t_dep_mean(7-i,:), 'o', 'LineWidth', 1);
hold on
end
text(0.05, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$T_{\mathrm{decrease}}$ [s]','Interpreter','Latex','fontsize',12);
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02', '$\tilde{\Theta}$=0.01', 'Interpreter','Latex', 'fontsize',10)

figure
for i = 1:5
errorbar(Shields(2:6), t_sat_mean(2:6,i),t_sat_std(2:6,i), 'o');
hold on
end
xlabel('$\tilde{\Theta}$','Interpreter','Latex','fontsize',12);ylabel('$T_{\mathrm{saturation}}$ [s]','Interpreter','Latex','fontsize',12);
legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);


%integrate u with t to get x vectors
tx = [0, data10.t(1,1:end)-1];
x=cell(6,5);Qx=cell(6,5);
for i = 1:6
    for j = 1:5
        x{i, j} = zeros(1, length(Q10)+1); 
        Qx{i, j} = zeros(1, length(Q10)+1); 
    end
end
for i=1:6
    for j=0:4
        eval(['Upx_int=[0, Upx',num2str(i),num2str(j),'];']);
        for k=2:length(Q10)+1
            eval(['x{',num2str(i),',',num2str(j+1),'}(',num2str(k),')=trapz(tx(1:',num2str(k),'), Upx_int(1:',num2str(k),'));']);
        end
    eval(['Qx{',num2str(i),',',num2str(j+1),'}=[Qi, fac_dim*Q',num2str(i),num2str(j),'];']);
    end
end

figure
for i=1:6
subplot(3,2,i);
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
box on
for j=0:4
   eval(['plot(x{',num2str(i),',',num2str(j+1),'}, Qx{',num2str(i),',',num2str(j+1),'});']); 
   hold on
end
ylim([0, 0.7]);
xlim([0, 20]);
xlabel('$x$ [m]','Interpreter','Latex','fontsize',12);ylabel('$Q$ [kg/m/s]','Interpreter','Latex','fontsize',12);
eval(['title("$\tilde{\Theta}$=0.0',num2str(i),'","Interpreter","Latex");']);
if i==1
   legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);
end
end

%calculate the peak point and relaxation length 
x_peak_mean=zeros(6,5);x_peak_std=zeros(6,5);x_sat_mean=zeros(6,5);x_sat_std=zeros(6,5);
for i=1:6
    for j=0:4
        eval(['[Qxpeak,x_peak, index_peak_x, x_peak_error]=getIndexPeak(x{',num2str(i),',',num2str(j+1),'}, Qx{',num2str(i),',',num2str(j+1),'});']);
        x_peak_mean(i,j+1)=x_peak;
        eval(['Qxpeak(',num2str(i),',',num2str(j+1),')=Qxpeak;']);
        %search from index_peak
        eval(['[x_sat,x_sat_error]=getSaturationIndex(x{',num2str(i),',',num2str(j+1),'}(index_peak_x:end),Qx{',num2str(i),',',num2str(j+1),'}(index_peak_x:end),Qs(',num2str(i),'));']);
        x_sat_mean(i,j+1)=x_sat-x_peak;
    end
end

figure
subplot(1,2,1)
for i = 1:6
plot(Liquid, x_peak_mean(7-i,:), 'o', 'LineWidth', 1);
hold on
end
ylim([0 18]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$L_{\mathrm{peak}}$ [m]','Interpreter','Latex','fontsize',12);
subplot(1,2,2)
for i = 1:6
plot(Liquid, x_sat_mean(7-i,:), 'o', 'LineWidth', 1);
hold on
end
text(0.05, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$L_{\mathrm{decrease}}$ [m]','Interpreter','Latex','fontsize',12);
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02', '$\tilde{\Theta}$=0.01','Interpreter','Latex', 'fontsize',10)


%linear fit to Lp=(A*omega+B)*(ustar-(C*omega+D))=(A*omega+B)*ustar-(A*omega+B)*(C*omega+D)
% a is (A*omega+B) and b is -(A*omega+B)*(C*omega+D); c is (C*omega+D)
a=zeros(1,5);b=zeros(1,5);
for i=1:size(x_peak_mean,2)
p=polyfit(ustar, x_peak_mean(:,i), 1);
a(i)=p(1);
b(i)=p(2);
end
c=-b./a;
Omega=Liquid*0.01;

p_a = polyfit(Omega, a, 1);
p_c = polyfit(Omega, c, 1);
b_fit = -polyval(p_a,Omega).*polyval(p_c,Omega);
Lp_fit=polyval(p_a,Omega)'*ustar+b_fit';
%van rijn and strypsteen (2020)
Fc=1.5e4*0.00025*(ustar-p_c(2))/(sqrt(D)*sqrt(9.81))+63*p_c(2);

%%plot the Lp-ustar, the fitted Lp-ustar with the inset including the
%%strypsteen data and formula, as well as the Lp-omega-Nip
colors_initial=[0 0 0;0 0 1;1 0 0];
figure
subplot1=subplot(1,2,1);
box on
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
for i = 1:size(x_peak_mean,2)
plot(ustar, x_peak_mean(:,i), 'o', 'LineWidth', 1);
hold on
end
for i = 1:size(Lp_fit,1)
plot(ustar, Lp_fit(i,:), '--', 'LineWidth', 1);
hold on
end
ylim([0 22]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$u_*$ [m/s]','Interpreter','Latex','fontsize',12);ylabel('$L_{\mathrm{peak}}$ [m]','Interpreter','Latex','fontsize',12);
legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);
subplot1_position = get(subplot1, 'Position');
inset_position = [subplot1_position(1) + 0.2 * subplot1_position(3), ... % left
                  subplot1_position(2) + 0.45 * subplot1_position(4), ... % bottom
                  0.35 * subplot1_position(3), ...                       % width
                  0.35 * subplot1_position(4)];                          % height
inset_axes = axes('Position', inset_position);
box(inset_axes, 'on');
% Plot data for the inset figure
set(gca, 'ColorOrder', [0 0 1;0 0 0;colors], 'NextPlot', 'replacechildren');
plot(inset_axes,[0.43, 0.475],[105, 82.5],'x', 'LineWidth', 1);
hold on;
plot(inset_axes, ustar,Fc,'-', 'LineWidth', 1); %'van rijn and strypsteen (2020)'
for i = 1:size(Lp_fit,1)
plot(inset_axes,ustar, Lp_fit(i,:), '--', 'LineWidth', 1);
hold on;
end
ylim(inset_axes, [0 120]);
xlabel(inset_axes, '$u_*$ [m/s]','Interpreter','Latex','fontsize',12);ylabel(inset_axes,'$L_{\mathrm{peak}}$ [m]','Interpreter','Latex','fontsize',12);
legend('Wet (Strypsteen et al., 2024)','Strypsteen et al. (2024) formula','This study');
%%plot Lpeak - Qi
subplot(1,2,2)
box on
set(gca, 'ColorOrder', colors_initial, 'NextPlot', 'replacechildren');
for i=1:3
    plot(Liquid, Lpeakthree(i,:),'o', 'LineWidth', 1);
    hold on;
end
ylim([0 55]);
text(0.05, 0.95, '(b)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);
ylabel('$L_\mathrm{peak}$ [m]','Interpreter','Latex','fontsize',12);
legend('$N_\mathrm{ip}$=1','$N_\mathrm{ip}$=5','$N_\mathrm{ip}$=10','Interpreter','Latex','fontsize',12);

function [t_steady, t_steady_error]=getSaturationIndex(t,array,ValueSteady)
array_smooth=array;
delta_y = array_smooth-ValueSteady;                % Compute differences between consecutive y values
tolerance = 0.15 * ValueSteady;        % Define a tolerance (1% of the final y value)

% Find the index where changes are smaller than the tolerance for a prolonged region
steady_idx = find(abs(delta_y) < tolerance, 1, 'first');  % First index meeting the criteria
t_steady = t(steady_idx);         % x-value where y reaches steady state
error_idx = find(abs(delta_y) < tolerance*5, 1, 'first');
% Error estimation
t_steady_error = (t_steady-t(error_idx))/2;
end