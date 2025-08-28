close all
clear
baseDir   = fullfile('stat','LBIniTP');
moistTags = {'Dry','M1','M5','M10','M20'};   % Omega order
Slist     = 1:6;                              % S001..S006 (Shields 0.01..0.06)
suffixes  = {'LBIni.stat'};    

data = cell(1, numel(Slist)*numel(moistTags));  % 1x30 cell

for s = 1:numel(Slist)                 % Shields index: 1..6
    S = Slist(s);
    for m = 1:numel(moistTags)         % Omega index: 1..5
        idx  = (s-1)*numel(moistTags) + m;  % linear position in data{...}
        stem = sprintf('S%03d%s', S, moistTags{m});

        % find existing file: try LBIni, then LBR1Ini
        path = '';
        for k = 1:numel(suffixes)
            cand = fullfile(baseDir, [stem suffixes{k}]);
            if exist(cand,'file'), path = cand; break; end
        end

        if isempty(path)
            warning('Missing file for %s (tried %s)', stem, strjoin(suffixes, ', '));
            data{idx} = [];
        else
            data{idx} = readMercuryCG(path);
        end
    end
end
data_ini=readMercuryCG('stat\LBIniTP\S005M20LBIniFirstTimeStep.stat');

D = 0.00025;
dz = 0.00025;%same as h in CG
Liquid=[0 1 5 10 20];
Shields=linspace(0.01,0.06,6);
ustar=sqrt(Shields*(2650-1.225)*9.81*D/1.225);
ustart=sqrt(0.005*(2650-1.225)*9.81*D/1.225);
dt=0.01;
fac_dim=sqrt((2650/1.225-1)*9.81*D^3)*2650;
colors = [0 0 0;0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0 0.4470 0.7410;0.6350 0.0780 0.1840];
num = 30; %5 Omega*6 Shields

Q    = cell(1, num);
U    = cell(1, num);
C    = cell(1, num);
H    = cell(1, num);
ids   = zeros(1, num);
for i=1:30
    id_shields = floor((i - 1)/5) + 1;
    [Q{i}, U{i}, C{i}, H{i}, Hsalsmooth] = CharTester(data{i}, id_shields, D, dz);
    ids(i) = steady_start_within15_afterpeak(C{i});
end
[Qi, Ui, Ci, Hi, Hisal] = CharTester(data_ini, 5, D, dz);
ids(14) = 250;
t = data{1}.t(1,:)-1;

figure;
for s = 1:6                         % 6 subplots for Shields 0.01..0.06
    subplot(3,2,s); hold on; grid on;
    for k = 1:5
        idx = (s-1)*5 + k;          % dataset index 1..30
        if idx > numel(C), break; end
        y = C{idx}(:).';            % ensure row
        plot(t, y, 'LineWidth', 1.2, 'Color', colors(k,:));

        % --- mark steady-start if valid ---
        if ~isempty(ids) && numel(ids) >= idx ...
                && ~isnan(ids(idx)) && ids(idx) >= 1 && ids(idx) <= numel(y)
            ti = t(ids(idx));
            yi = y(ids(idx));
            plot(ti, yi, 'o', ...
                'MarkerSize', 6, ...
                'MarkerEdgeColor', colors(k,:), ...
                'MarkerFaceColor', colors(k,:));
            % Optional label:
            % text(ti, yi, sprintf('  %d', ids(idx)), 'Color', colors(k,:), 'VerticalAlignment','middle');
            % Optional vertical guide:
            % xline(ti, '--', 'Color', colors(k,:), 'HandleVisibility','off');
        end
    end
    title(sprintf('\\Theta = %.2f', Shields(s)));
    xlabel('t'); ylabel('C');
end

Qs = nan(6,5);
Qstd = nan(6,5);
Cs = nan(6,5);
Us = nan(6,5);
Hs = nan(6,5);
for i = 1:num
    s = ceil(i / 5);            % Shields row (1..6)
    o = mod(i-1, 5) + 1;        % Omega  col (1..5)

    % guard invalid ids
    if isempty(ids) || isnan(ids(i)) || ids(i) < 1 || ids(i) > numel(Q{i})
        continue;                 % leave as NaN
    end

    idx0 = ids(i):numel(Q{i});    % steady segment i..end

    Qs(s,o) = mean(Q{i}(idx0));
    Qstd(s,o) = std(Q{i}(idx0));
    Cs(s,o) = mean(C{i}(idx0));
    Us(s,o) = mean(U{i}(idx0));
    Hs(s,o) = mean(H{i}(idx0));
end

figure
subplot(2,2,1)
for i=1:6
    plot(Liquid, Qs(7-i,:),'o','LineWidth', 1);
    hold on
end
plot(Liquid,mean(Qs),'k-');
ylim([0 0.05]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$Q_{\mathrm{steady}}$ [kg/m/s]','Interpreter','Latex','fontsize',12);
subplot(2,2,2)
for i=1:6
    plot(Liquid, Cs(7-i,:),'o','LineWidth', 1);
    hold on
end
plot(Liquid,mean(Cs),'k-');
ylim([0 0.06]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$C_{\mathrm{steady}}$ [kg/m$^2$]','Interpreter','Latex','fontsize',12);
subplot(2,2,3)
for i=1:6
    plot(Liquid, Us(7-i,:),'o','LineWidth', 1);
    hold on
end
plot(Liquid,mean(Us),'k-');
ylim([0 3.5]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$U_{\mathrm{steady}}$ [m/s]','Interpreter','Latex','fontsize',12);
subplot(2,2,4)
for i=1:6
    plot(Liquid, Hs(7-i,:),'o','LineWidth', 1);
    hold on
end
plot(Liquid,mean(Hs),'k-');
ylim([0 0.03]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
xlabel('$\Omega$ [$\%$]','Interpreter','Latex','fontsize',12);ylabel('$H_{\mathrm{steady}}$ [m]','Interpreter','Latex','fontsize',12);
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02', '$\tilde{\Theta}$=0.01', 'Interpreter','Latex', 'fontsize',10);

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
p = polyfit(Shields, Qs(:,1), 1); 
x_fit = linspace(min(Shields), max(Shields), 100); 
y_fit = polyval(p, x_fit); 

%%validation figure for dry limit
%dry Q-t
figure
subplot(1,2,1)
for i=1:6
    id = 26 - (i-1)*5;  
    plot([0, t], [Qi(1), fac_dim*Q{id}]);
    hold on
end
ylim([0 0.7]);
text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 12);
legend('$\tilde{\Theta}$=0.06','$\tilde{\Theta}$=0.05','$\tilde{\Theta}$=0.04','$\tilde{\Theta}$=0.03','$\tilde{\Theta}$=0.02','$\tilde{\Theta}$=0.01','Interpreter','Latex', 'fontsize',12)
xlabel('$t$ [s]','Interpreter','Latex','fontsize',12);ylabel('$Q$ [kg/m/s]','Interpreter','Latex','fontsize',12);
set(gca,'FontSize',14);
subplot(1,2,2)
errorbar(Shields, Qs(:,1), Qstd(:,1), 'or');
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
ylim([-0.0001 0.045]);


function idx0 = steady_start_within15_afterpeak(x)
% Return the starting index (1-based) of the last-longest segment whose
% values stay within 10% of the average of the last 100 samples of x,
% constrained to start after the global peak.
%
%   idx0 = steady_start_within10_afterpeak(x)

    x = x(:);
    N = numel(x);
    if N == 0
        idx0 = NaN; return;
    end

    % --- average over the last up-to-100 samples ---
    k = min(150, N);
    xbar = mean(x(N-k+1:N));

    % --- 15% relative band around that average ---
    tol = 0.15 * abs(xbar);
    if tol == 0
        mask = (x == 0);
    else
        mask = abs(x - xbar) <= tol;
    end

    % --- enforce "after the peak" (use last occurrence of the global max) ---
    [~, imax_first] = max(x);
    imax = find(x == x(imax_first), 1, 'last');  % robust to flat maxima
    mask(1:imax) = false;

    % --- find contiguous true-runs in mask ---
    d = diff([false; mask; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    if isempty(starts)
        idx0 = NaN;  % no qualifying segment after the peak
        return;
    end

    % --- choose the last-longest segment (tie-break favors later start) ---
    lens = ends - starts + 1;
    [~, kbest] = max(lens + starts/1e9);
    idx0 = starts(kbest);
end

