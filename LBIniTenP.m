close all
% clear
% data10=readMercuryCG('stat\LBIniTP\S001DryLBIni.stat');
% data11=readMercuryCG('stat\LBIniTP\S001M1LBIni.stat');
% data12=readMercuryCG('stat\LBIniTP\S001M5LBIni.stat');
% data13=readMercuryCG('stat\LBIniTP\S001M10LBIni.stat');
% data14=readMercuryCG('stat\LBIniTP\S001M20LBIni.stat');
% 
% data20=readMercuryCG('stat\LBIniTP\S002DryLBIni.stat');
% data21=readMercuryCG('stat\LBIniTP\S002M1LBIni.stat');
% data22=readMercuryCG('stat\LBIniTP\S002M5LBIni.stat');
% data23=readMercuryCG('stat\LBIniTP\S002M10LBIni.stat');
% data24=readMercuryCG('stat\LBIniTP\S002M20LBIni.stat');
% 
% data30=readMercuryCG('stat\LBIniTP\S003DryLBIni.stat');
% data31=readMercuryCG('stat\LBIniTP\S003M1LBIni.stat');
% data32=readMercuryCG('stat\LBIniTP\S003M5LBIni.stat');
% data33=readMercuryCG('stat\LBIniTP\S003M10LBIni.stat');
% data34=readMercuryCG('stat\LBIniTP\S003M20LBIni.stat');
%  
% data40=readMercuryCG('stat\LBIniTP\S004DryLBIni.stat');
% data41=readMercuryCG('stat\LBIniTP\S004M1LBR1Ini.stat');
% data42=readMercuryCG('stat\LBIniTP\S004M5LBIni.stat');
% data43=readMercuryCG('stat\LBIniTP\S004M10LBIni.stat');
% data44=readMercuryCG('stat\LBIniTP\S004M20LBIni.stat');
% 
% data50=readMercuryCG('stat\LBIniTP\S005DryLBIni.stat');
% data51=readMercuryCG('stat\LBIniTP\S005M1LBIni.stat');
% data52=readMercuryCG('stat\LBIniTP\S005M5LBR1Ini.stat');
% data53=readMercuryCG('stat\LBIniTP\S005M10LBR1Ini.stat');
% data54=readMercuryCG('stat\LBIniTP\S005M20LBIni.stat');
% 
% data60=readMercuryCG('stat\LBIniTP\S006DryLBIni.stat');
% data61=readMercuryCG('stat\LBIniTP\S006M1LBIni.stat');
% data62=readMercuryCG('stat\LBIniTP\S006M5LBIni.stat');
% data63=readMercuryCG('stat\LBIniTP\S006M10LBIni.stat');
% data64=readMercuryCG('stat\LBIniTP\S006M20LBIni.stat');
% data_ini=readMercuryCG('stat\LBIniTP\S005M20LBIniFirstTimeStep.stat');

D = 0.00025;
dz = 0.00025;%same as h in CG
Liquid=[0 1 5 10 20];
Shields=linspace(0.01,0.06,6);
%for validation
Dexp=0.000175;
Ufs=linspace(10,16,4);
u_star_exp=0.4*Ufs/log(0.6/(Dexp/60));%assume z0=d/60
Shields_exp=u_star_exp.^2*1.225/(2650-1.225)/9.81/Dexp;%density is 2650, in paper 1520 bulk
Liquid_exp=[0.143 0.351 1.443 2.713]*2.65;
Q_exp=[0.2 0.25 0.45 0.2; 0.7 0.75 0.6 0.21; 1.3 1.4 0.7 0.22; 2.2 2.6 0.8 0.25]*0.001*100/sqrt((2650/1.225-1)*9.81*Dexp^3)/2650;

for j=1:6
for i=0:4
    eval(['[Q',num2str(j),num2str(i),', Upx',num2str(j), num2str(i), ', Cvx',num2str(j), num2str(i),', H',num2str(j), num2str(i),'] = CharTester(data',num2str(j), num2str(i),',',num2str(j),',D, dz);']); 
end
end
[Qini,Upini,Cvxini,Hini]=CharTester(data_ini,5,D,dz);
Cini=Cvxini(1);

%Qs,Qstd,Cs,Ups
fac_dim=sqrt((2650/1.225-1)*9.81*D^3)*2650;
for i=1:6
for j=0:4
if i<4
   eval(['Qs(',num2str(i),',',num2str(j+1),')=fac_dim*mean(Q',num2str(i),num2str(j),'(2.5/5*end:end));']);
   eval(['Cs(',num2str(i),',',num2str(j+1),')=mean(Cvx',num2str(i),num2str(j),'(2.5/5*end:end));']);
   eval(['Ups(',num2str(i),',',num2str(j+1),')=getMeanOfNonZero(Upx',num2str(i),num2str(j),'(2.5/5*end:end));']);
   eval(['Hs(',num2str(i),',',num2str(j+1),')=getMeanOfNonZero(H',num2str(i),num2str(j),'(2.5/5*end:end));']);
   else
       eval(['Qs(',num2str(i),',',num2str(j+1),')=fac_dim*mean(Q',num2str(i),num2str(j),'(3.5/5*end:end));']); 
       eval(['Cs(',num2str(i),',',num2str(j+1),')=mean(Cvx',num2str(i),num2str(j),'(3.5/5*end:end));']);
       eval(['Ups(',num2str(i),',',num2str(j+1),')=getMeanOfNonZero(Upx',num2str(i),num2str(j),'(3.5/5*end:end));']);
       eval(['Hs(',num2str(i),',',num2str(j+1),')=getMeanOfNonZero(H',num2str(i),num2str(j),'(3.5/5*end:end));']);
end
end
end

%%define the transient phase by C-t
for i=1:6
    for j=0:4
        eval(['index_peak=getIndexPeak(Cvx',num2str(i),num2str(j),');']);
        t_peak(i,j+1)=index_peak*0.01;
        %search from index_peak
        eval(['index_tran=getTtranIndex(Cvx',num2str(i),num2str(j),'(index_peak:end),Cs(',num2str(i),',',num2str(j+1),'));']);
        t_tran(i,j+1)=(index_peak+index_tran-1)*0.01;
    end
end
%t_tran(1,2)=2.67;

for i=1:6
    for j=0:4
        eval(['Qero(',num2str(i),',',num2str(j+1),')=fac_dim*mean(Q',num2str(i),num2str(j),'(1:t_peak(',num2str(i),',',num2str(j+1),')*100));']); 
        eval(['Cero(',num2str(i),',',num2str(j+1),')=mean(Cvx',num2str(i),num2str(j),'(1:t_peak(',num2str(i),',',num2str(j+1),')*100));']); 
        eval(['Upero(',num2str(i),',',num2str(j+1),')=getMeanOfNonZero(Upx',num2str(i),num2str(j),'(1:t_peak(',num2str(i),',',num2str(j+1),')*100));']); 
        eval(['Hero(',num2str(i),',',num2str(j+1),')=getMeanOfNonZero(H',num2str(i),num2str(j),'(1:t_peak(',num2str(i),',',num2str(j+1),')*100));']); 
        
        eval(['Qdep(',num2str(i),',',num2str(j+1),')=fac_dim*mean(Q',num2str(i),num2str(j),'(t_peak(',num2str(i),',',num2str(j+1),')*100:t_tran(',num2str(i),',',num2str(j+1),')*100));']); 
        eval(['Cdep(',num2str(i),',',num2str(j+1),')=mean(Cvx',num2str(i),num2str(j),'(t_peak(',num2str(i),',',num2str(j+1),')*100:t_tran(',num2str(i),',',num2str(j+1),')*100));']); 
        eval(['Updep(',num2str(i),',',num2str(j+1),')=mean(Upx',num2str(i),num2str(j),'(t_peak(',num2str(i),',',num2str(j+1),')*100:t_tran(',num2str(i),',',num2str(j+1),')*100));']); 
        eval(['Hdep(',num2str(i),',',num2str(j+1),')=mean(H',num2str(i),num2str(j),'(t_peak(',num2str(i),',',num2str(j+1),')*100:t_tran(',num2str(i),',',num2str(j+1),')*100));']); 
    end
end

%%%%%%plotting
%% time series
figure
for i=1:6
subplot(3,2,i);
eval(['plot(data',num2str(i),'0.t(1,:)-1, fac_dim*Q',num2str(i),'0);']); 
hold on
for j=1:4
   eval(['plot(data',num2str(i),num2str(j),'.t(1,:)-2, fac_dim*Q',num2str(i),num2str(j),');']); 
   hold on
end
ylim([0, 0.7]);xlim([0, 5]);
xlabel('$t$ [s]','Interpreter','Latex','fontsize',12);ylabel('$Q$ [kg/m/s]','Interpreter','Latex','fontsize',12);
eval(['title("$\tilde{\Theta}$=0.0',num2str(i),'","Interpreter","Latex");']);
if i==1
   legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);
end
end

% figure
% for i=1:6
% subplot(3,2,i);
% eval(['plot([0 data',num2str(i),'0.t(1,:)-1], [Cini Cvx',num2str(i),'0]);']); 
% hold on
% for j=1:4
%    eval(['plot([0 data',num2str(i),num2str(j),'.t(1,:)-2], [Cini Cvx',num2str(i),num2str(j),']);']); 
%    hold on
% end
% ylim([0, 0.6]);xlim([0, 5]);
% xlabel('$t$ [s]','Interpreter','Latex','fontsize',12);ylabel('$C_\mathrm{sal}$ [kg/m$^2$]','Interpreter','Latex','fontsize',12);
% eval(['title("$\tilde{\Theta}$=0.0',num2str(i),'","Interpreter","Latex");']);
% if i==1
%    legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);
% end
% end

% figure
% for i=1:6
% subplot(3,2,i);
% eval(['plot(data',num2str(i),'0.t(1,:)-1, Upx',num2str(i),'0);']); 
% hold on
% for j=1:4
%    eval(['plot(data',num2str(i),num2str(j),'.t(1,:)-2, Upx',num2str(i),num2str(j),');']); 
%    hold on
% end
% ylim([0, 10]);xlim([0, 5]);
% xlabel('$t$ [s]','Interpreter','Latex','fontsize',12);ylabel('$U_\mathrm{sal}$ [m/s]','Interpreter','Latex','fontsize',12);
% eval(['title("$\tilde{\Theta}$=0.0',num2str(i),'","Interpreter","Latex");']);
% if i==1
%    legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$','Interpreter','Latex', 'fontsize',10);
% end
% end

%% eroding and depositing time scale
% figure
% subplot(1,2,1)
% for i = 1:length(Shields)
% plot(Liquid, t_peak(length(Shields)-i+1,:), 'o--', 'LineWidth', 1.5, 'MarkerSize', 8);
% hold on
% end
% ylim([0 3.5]);
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$t_{\mathrm{eroding}}$ [s]','Interpreter','Latex','fontsize',12);
% subplot(1,2,2)
% for i = 1:length(Shields)
% plot(Liquid, t_tran(length(Shields)-i+1,:)-t_peak(length(Shields)-i+1,:), 'o--', 'LineWidth', 1.5, 'MarkerSize', 8);
% hold on
% end
% ylim([0 3.5]);
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$t_{\mathrm{depositing}}$ [s]','Interpreter','Latex','fontsize',12);
% legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02','$\tilde{\Theta}$=0.01', 'Interpreter','Latex', 'fontsize',10)
% 

% figure
% subplot(1,4,1)
% for i=1:6
% plot(Liquid,Qero(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
% hold on
% end
% title('Simulation results');
% ylim([0 0.4]);
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\tilde{Q}_{\mathrm{eroding}}$ [-]','Interpreter','Latex','fontsize',12);
% legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02','$\tilde{\Theta}$=0.01', 'Interpreter','Latex', 'fontsize',10)
% subplot(1,4,2)
% for i=1:6
% plot(Liquid,Qdep(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
% hold on
% end
% title('Simulation results');
% ylim([0 0.4]);
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\tilde{Q}_{\mathrm{depositing}}$ [-]','Interpreter','Latex','fontsize',12);
% subplot(1,4,3)
% for i=1:6
% plot(Liquid,Qs(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
% hold on
% end
% title('Simulation results');
% ylim([0 0.4]);
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\tilde{Q}_{\mathrm{steady}}$ [-]','Interpreter','Latex','fontsize',12);


figure
subplot(3,3,1)
for i=1:6
plot(Liquid,Qero(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 0.45]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$Q_{\mathrm{eroding}}$ [kg/m/s]','Interpreter','Latex','fontsize',12);
title('Eroding phase');
subplot(3,3,2)
for i=1:6
plot(Liquid,Qdep(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 0.45]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$Q_{\mathrm{depositing}}$ [kg/m/s]','Interpreter','Latex','fontsize',12);
title('Depositing phase');
subplot(3,3,3)
for i=1:6
plot(Liquid,Qs(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 0.45]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$Q_{\mathrm{steady}}$ [kg/m/s]','Interpreter','Latex','fontsize',12);
title('Steady state');
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02','$\tilde{\Theta}$=0.01', 'Interpreter','Latex', 'fontsize',10)
subplot(3,3,4)
for i=1:6
plot(Liquid,Cero(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 0.3]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$C_{\mathrm{sal,eroding}}$ [kg/m$^2$]','Interpreter','Latex','fontsize',12);
subplot(3,3,5)
for i=1:6
plot(Liquid,Cdep(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 0.3]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$C_{\mathrm{sal,depositing}}$ [kg/m$^2$]','Interpreter','Latex','fontsize',12);
subplot(3,3,6)
for i=1:6
plot(Liquid,Cs(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 0.3]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$C_{\mathrm{sal,steady}}$ [kg/m$^2$]','Interpreter','Latex','fontsize',12);
subplot(3,3,7)
for i=1:6
plot(Liquid,Upero(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 6.1]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$U_{\mathrm{sal,eroding}}$ [m/s]','Interpreter','Latex','fontsize',12);
subplot(3,3,8)
for i=1:6
plot(Liquid,Updep(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 6.1]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$U_{\mathrm{sal,depositing}}$ [m/s]','Interpreter','Latex','fontsize',12);
subplot(3,3,9)
for i=1:6
plot(Liquid,Ups(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
hold on
end
ylim([0 6.1]);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$U_{\mathrm{sal,steady}}$ [m/s]','Interpreter','Latex','fontsize',12);





%%compare with Han et al. (2011)
%get mean values of all cases by the same time scale (tero(6,1))
% for i=1:6
%     for j=0:4
%         eval(['Qsim(',num2str(i),',',num2str(j+1),')=mean(Q',num2str(i),num2str(j),'(1:t_peak(6,1)*100));']); 
%     end
% end
% 
% figure
% subplot(1,2,1)
% for i=1:6
% plot(Liquid,Qsim(length(Shields)-i+1,:),'--o','LineWidth', 1.5, 'MarkerSize', 8);
% hold on
% end
% ylim([0 0.4]);
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\tilde{Q}_{\mathrm{simulation}}$ [-]','Interpreter','Latex','fontsize',12);
% legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02','$\tilde{\Theta}$=0.01', 'Interpreter','Latex', 'fontsize',10)
% subplot(1,2,2)
% for i=1:size(Q_exp,1)
% plot(Liquid_exp, Q_exp(size(Q_exp,1)-i+1,:),'--d','LineWidth', 1.5, 'MarkerSize', 8);
% hold on
% end
% ylim([0 0.4]);box on
% xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$\tilde{Q}_{\mathrm{experiment}}$ [-]','Interpreter','Latex','fontsize',12);
% legend('$\tilde{\Theta}$=0.0737','$\tilde{\Theta}$=0.0564','$\tilde{\Theta}$=0.0415','$\tilde{\Theta}$=0.0288', 'Interpreter','Latex', 'fontsize',10)



function index_peak=getIndexPeak(Value)
Value_smooth=smoothdata(Value);
[max_value, max_index] = max(Value_smooth);
index_peak=max_index;
end

function index_tran=getTtranIndex(array,ValueSteady)
array_smooth=smoothdata(array);
value = ValueSteady;
threshold = 0.05*value;
if threshold <= 1e-4
    index_tran = find(array_smooth - value <= threshold, 1, 'first');
else
    index_tran = find(abs(array_smooth - value) <= threshold, 1, 'first');
end
end

function smoothed_data=smoothdata(Q)
data = Q;
% Define window size for moving average
window_size = 30;
% Apply moving average filter
smoothed_data = movmean(data, window_size);
end
