close all
% clear
% data0=readMercuryCG('stat\LBIniTP\zprofile\S005DryLBIniW50Dh5D.stat');
% data1=readMercuryCG('stat\LBIniTP\zprofile\S005M1LBIniW50Dh5D.stat');
% data2=readMercuryCG('stat\LBIniTP\zprofile\S005M5LBR1IniW50Dh5D.stat');
% data3=readMercuryCG('stat\LBIniTP\zprofile\S005M10LBR1IniW50Dh5D.stat');
% data4=readMercuryCG('stat\LBIniTP\zprofile\S005M20LBIniW50Dh5D.stat');
% 
% data5=readMercuryCG('stat\LBIniTP\zprofile\S006DryLBIniW50Dh5D.stat');
% data6=readMercuryCG('stat\LBIniTP\zprofile\S006M1LBIniW50Dh5D.stat');
% data7=readMercuryCG('stat\LBIniTP\zprofile\S006M5LBIniW50Dh5D.stat');
% data8=readMercuryCG('stat\LBIniTP\zprofile\S006M10LBIniW50Dh5D.stat');
% data9=readMercuryCG('stat\LBIniTP\zprofile\S006M20LBIniW50Dh5D.stat');

% data10=readMercuryCG('stat\LBIniTP\zprofile\S004DryLBIniW50Dh5D.stat');
% data11=readMercuryCG('stat\LBIniTP\zprofile\S004M1LBR1IniW50Dh5D.stat');
% data12=readMercuryCG('stat\LBIniTP\zprofile\S004M5LBIniW50Dh5D.stat');
% data13=readMercuryCG('stat\LBIniTP\zprofile\S004M10LBIniW50Dh5D.stat');
% data14=readMercuryCG('stat\LBIniTP\zprofile\S004M20LBIniW50Dh5D.stat');

D=0.00025;
dz=5*D;
z=linspace(0, 0.2, 160);
Liquid=[0 1 5 10 20];
t_ero=[0.69 0.85 1.13 1.48 2.65 0.62 0.79 1.53 1.91 2.44]*100;
t_tran=t_ero+[1.05 1.1 0.79 0.68 0.54 1.05 0.96 0.6 1.73 0.98]*100;

%%vertical flux profile
for i=0:9
eval(['[Cz',num2str(i),',Upz',num2str(i),',Qz',num2str(i),'] = getZProfile(data',num2str(i),', dz);']); 
eval(['Qzero(',num2str(i+1),',:)=mean(Qz',num2str(i),'(:,1:t_ero(',num2str(i+1),')),2);']);
eval(['Czero(',num2str(i+1),',:)=mean(Cz',num2str(i),'(:,1:t_ero(',num2str(i+1),')),2);']);
eval(['Upzero(',num2str(i+1),',:)=getMeanOfNonZero(Upz',num2str(i),'(:,1:t_ero(',num2str(i+1),')));']);
eval(['Czs(',num2str(i+1),',:)=mean(Cz',num2str(i),'(:,400:end),2);']);
eval(['Upzs(',num2str(i+1),',:)=getMeanOfNonZero(Upz',num2str(i),'(:,400:end));']);
eval(['Qzs(',num2str(i+1),',:)=mean(Qz',num2str(i),'(:,400:end),2);']);
end

y_grey=12*D+50*D;
figure
subplot(1,3,1)
for i=1:5
plot(Qzero(i,:), z);
hold on
end
% x=[0 0.02 0.02 0];y=[0 0 y_grey y_grey];
% fill(x,y,[.7 .7 .7],'FaceAlpha', 0.5);
ylim([y_grey 0.1]);
xlabel('$Q$(z) [kg/m/s]','Interpreter','Latex','fontsize',12);ylabel('$z$ [m]','Interpreter','Latex','fontsize',12);
subplot(1,3,2)
for i=1:5
semilogx(Czero(i,:), z);
hold on
end
% x=[1e-5 1e0 1e0 1e-5];y=[0 0 y_grey y_grey];
% fill(x,y,[.7 .7 .7],'FaceAlpha', 0.5);
ylim([y_grey 0.1]);
xlabel('$C_\mathrm{sal}$(z) [kg/m$^2$]','Interpreter','Latex','fontsize',12);ylabel('$z$ [m]','Interpreter','Latex','fontsize',12);
%title('$\tilde{\Theta}$=0.05, Eroding phase','Interpreter','Latex', 'fontsize',12);
subplot(1,3,3)
for i=1:5
plot(Upzero(i,:), z);
hold on
end
% x=[0 10 10 0];y=[0 0 y_grey y_grey];
% fill(x,y,[.7 .7 .7],'FaceAlpha', 0.5);
xlim([0 10]);ylim([y_grey 0.1]);
xlabel('$U_\mathrm{sal}$(z) [m/s]','Interpreter','Latex','fontsize',12);ylabel('$z$ [m]','Interpreter','Latex','fontsize',12);
legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$', 'Interpreter','Latex', 'fontsize',10);

figure
subplot(1,3,1)
for i=6:10
plot(Qzero(i,:), z);
hold on
end
ylim([y_grey 0.1]);
xlabel('$Q$(z) [kg/m/s]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
subplot(1,3,2)
for i=6:10
semilogx(Czero(i,:), z);
hold on
end
ylim([y_grey 0.1]);
xlabel('$C_\mathrm{sal}$(z) [kg/m$^2$]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
title('$\tilde{\Theta}$=0.06, Eroding phase','Interpreter','Latex', 'fontsize',12);
subplot(1,3,3)
for i=6:10
plot(Upzero(i,:), z);
hold on
end
xlim([0 12]);ylim([y_grey 0.1]);
xlabel('$U_\mathrm{sal}$(z) [m/s]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$', 'Interpreter','Latex', 'fontsize',10);

% figure
% subplot(1,3,1)
% for i=10:14
% plot(Qzero(i,:), z);
% hold on
% end
% x=[0 0.03 0.03 0];y=[0 0 y_grey y_grey];
% fill(x,y,[.7 .7 .7],'FaceAlpha', 0.5);
% % hold on
% % x_range=linspace(0, 0.03);
% % plot(x_range, 0.045*ones(1,length(x_range)),'--k');
% ylim([0 0.1]);
% xlabel('$Q$(z) [kg/m/s]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
% subplot(1,3,2)
% for i=10:14
% semilogx(Czero(i,:), z);
% hold on
% end
% x=[5e-7 1.09e0 1.09e0 5e-7];y=[0 0 y_grey y_grey];
% fill(x,y,[.7 .7 .7],'FaceAlpha', 0.5);
% % hold on
% % x_range=linspace(5e-7, 1.09e0);
% % plot(x_range, 0.045*ones(1,length(x_range)),'--k');
% ylim([0 0.1]);
% xlabel('$C_\mathrm{sal}$(z) [kg/m$^2$]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
% title('$\tilde{\Theta}$=0.04, Eroding phase','Interpreter','Latex', 'fontsize',12);
% subplot(1,3,3)
% for i=10:14
% plot(Upzero(i,:), z);
% hold on
% end
% x=[0 12 12 0];y=[0 0 y_grey y_grey];
% fill(x,y,[.7 .7 .7],'FaceAlpha', 0.5);
% % hold on
% % x_range=linspace(0,12);
% % plot(x_range, 0.045*ones(1,length(x_range)),'--k');
% xlim([0 12]);ylim([0 0.1]);
% xlabel('$U_\mathrm{sal}$(z) [m/s]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
% legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$', 'Interpreter','Latex', 'fontsize',10);
% 


figure
subplot(1,3,1)
for i=1:5
plot(Qzs(i,:), z);
hold on
end
ylim([y_grey 0.1]);
xlabel('$Q$(z) [kg/m/s]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
subplot(1,3,2)
for i=1:5
semilogx(Czs(i,:), z);
hold on
end
ylim([y_grey 0.1]);
xlabel('$C_\mathrm{sal}$(z) [kg/m$^2$]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
title('$\tilde{\Theta}$=0.05, Steady state','Interpreter','Latex', 'fontsize',12);
subplot(1,3,3)
for i=1:5
plot(Upzs(i,:), z);
hold on
end
ylim([y_grey 0.1]);
xlabel('$U_\mathrm{sal}$(z) [m/s]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$', 'Interpreter','Latex', 'fontsize',10);

figure
subplot(1,3,1)
for i=6:10
plot(Qzs(i,:), z);
hold on
end
ylim([y_grey 0.1]);
xlabel('$Q$(z) [kg/m/s]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
subplot(1,3,2)
for i=6:10
semilogx(Czs(i,:), z);
hold on
end
ylim([y_grey 0.1]);
xlabel('$C_\mathrm{sal}$(z) [kg/m$^2$]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
title('$\tilde{\Theta}$=0.06, Steady state','Interpreter','Latex', 'fontsize',12);
subplot(1,3,3)
for i=6:10
plot(Upzs(i,:), z);
hold on
end
ylim([y_grey 0.1]);
xlabel('$U_\mathrm{sal}$(z) [m/s]','Interpreter','Latex','fontsize',12);ylabel('z [m]','fontsize',12);
legend('$\Omega$=0','$\Omega=1\%$','$\Omega=5\%$','$\Omega=10\%$','$\Omega=20\%$', 'Interpreter','Latex', 'fontsize',10);


% figure
% for i=1:5:6
% plot(Upzero(i,:), z);
% hold on
% end
% xlabel('Upx(z) [m/s]');ylabel('$z/d_{mean}$ [-]','Interpreter','Latex','fontsize',12);
% legend('$\tilde{\Theta}$=0.05, $\Omega$=0','$\tilde{\Theta}$=0.06, $\Omega$=0','Interpreter','Latex', 'fontsize',10)
% figure
% for i=2:5:7
% plot(Upzero(i,:), z);
% hold on
% end
% legend('$\tilde{\Theta}$=0.05, $\Omega$=1$\%$','$\tilde{\Theta}$=0.06, $\Omega$=1$\%$', 'Interpreter','Latex', 'fontsize',10)
% figure
% for i=3:5:8
% plot(Upzero(i,:), z);
% hold on
% end
% legend('$\tilde{\Theta}$=0.05, $\Omega$=5$\%$','$\tilde{\Theta}$=0.06, $\Omega$=5$\%$', 'Interpreter','Latex', 'fontsize',10)
% figure
% for i=4:5:9
% plot(Upzero(i,:), z);
% hold on
% end
% legend('$\tilde{\Theta}$=0.05, $\Omega$=10$\%$','$\tilde{\Theta}$=0.06, $\Omega$=10$\%$', 'Interpreter','Latex', 'fontsize',10)



function Value_in=GetMeanInIncrease(Value,IsU)
Value_in=[];
for i=1:size(Value,1)
Value_smooth=smoothdata(Value(i,:));
[max_value, max_index] = max(Value_smooth);
if IsU == 1
   Value_in(i)=getMeanOfNonZero(Value(i,1:max_index));
else
   Value_in(i)=mean(Value(i,1:max_index));
end
end
end

function smoothed_data=smoothdata(Q)
data = Q;
% Define window size for moving average
window_size = 20;
% Apply moving average filter
smoothed_data = movmean(data, window_size);
end