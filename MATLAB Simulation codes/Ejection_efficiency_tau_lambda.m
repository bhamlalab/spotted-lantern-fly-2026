%% =================== Two-spring sweep + tau & Lambda + energy transfer + PRA alpha ===================
clear; close all; clc

% ---------- Base parameters ----------
pars.mdroplet = 1;
pars.kdroplet = 1;
zeta = 0;
pars.c = 2*zeta;

% Paper's definition: fo = (1/pi)*sqrt(kd/md)   (Rayleigh 2nd mode)
pars.fo = sqrt(pars.kdroplet/pars.mdroplet) / pi;
% 만약 기존 코드처럼 1/(2*pi)를 원하면 위 줄 대신:
 %pars.fo = sqrt(pars.kdroplet/pars.mdroplet) / (2*pi);

% Initial conditions
x0stylus  = -1;
x0droplet = -1;

% Data holders
F = []; Alpha = []; Lambda = [];
Tmax_vel = []; Tmax_comp = [];
FoF = []; Tau = []; fmin_vec = [];

% (추가) Energy transfer 저장소
Wraw = [];

% ---------- Sweep over i (changes stylus stiffness) ----------
for i = 1.0e-4 : 2.5e-3 : 10
    % Stylus params (샤프슈터 가정)
    pars.mstylus = 2000*pars.mdroplet;
    pars.kstylus = 2000*pars.kdroplet/(i);

    % Stylus driving frequency (engine)
    pars.fstylus = sqrt(pars.kstylus/(pars.mdroplet+pars.mstylus)) / (2*pi);

    % 대표 주파수(시뮬레이션 시간축 설정용): min(fo, fstylus)
    fmin = min(pars.fo, pars.fstylus);
    fmin_vec = [fmin_vec fmin];

    Tmin = 1/fmin;
    numcycles = 0.5;                          % integrate for half a min-period
    time = linspace(0, numcycles*Tmin, 30000);

    % ---------- ODE simulation ----------
    BC = [x0stylus 0 x0droplet 0];
    [t, dyn] = ode45(@supertwospring, time, BC, [], pars);

    xstylus  = dyn(:,1);  vstylus  = dyn(:,2);
    xdroplet = dyn(:,3);  vdroplet = dyn(:,4);

    x_rel = xdroplet - xstylus;               % relative compression (drop - stylus)

    % ---------- Extract timings ----------
    % tv: time of max stylus velocity
    [~, idx_vmax_s] = max(vstylus);
    tv = t(idx_vmax_s);

    % tc: time of first compression peak (max of -x_rel)
    [~, loc_comp] = findpeaks(-x_rel);
    if isempty(loc_comp)
        % 압축 피크를 못 찾으면 이번 포인트는 건너뜀
        continue
    end
    tc = t(loc_comp(1));

    % ---------- Metrics ----------
    % Lambda = Vd,max / Vs,max
    lambda = max(vdroplet) / max(vstylus);
    alpha  = lambda^2;

    % Tau (논문 정의): (tc - tv) / T  with  T = 1/fstylus
    tau = (tc - tv) * pars.fstylus;           % 부호 유지, 무차원

    % ---------- Store ----------
    F         = [F pars.fstylus];
    Alpha     = [Alpha alpha];
    Lambda    = [Lambda lambda];
    Tmax_vel  = [Tmax_vel tv];
    Tmax_comp = [Tmax_comp tc];
    FoF       = [FoF pars.fo/pars.fstylus];
    Tau       = [Tau tau];

    % ---------- (추가) Energy transfer (work) ----------
    % te: droplet(상부 질량) 속도가 최대가 되는 시점
    [~, idx_vmax_d] = max(vdroplet);
    te = t(idx_vmax_d);

    % Fd = kd * (xd - xs) : lower spring이 upper mass에 가하는 힘
    Fd_vec = pars.kdroplet * (xdroplet - xstylus);

    % W = \int_0^{te} Fd(t) * Vs(t) dt
    mask_te = t <= te;
    W_te = trapz(t(mask_te), Fd_vec(mask_te) .* vstylus(mask_te));

    % Store raw work
    Wraw = [Wraw W_te];
end

% ---------- (추가) PR Applied 정의의 ejection efficiency alpha ----------
% r = f_o/f
r = FoF(:);

% Ve/Vp* = r/(r-1) * sin( 2*pi/(r+1) )
Ve_over_Vp = r./(r-2) .* sin( 4*pi./(r+2) );

% r ~= 1 특이점 처리
Ve_over_Vp(abs(r-1) < 1e-6) = NaN;

% alpha_PRA = (Ve/Vp*)^2
Alpha_PRA = (Ve_over_Vp).^2;

% =================== Plots ===================

% Fig A — tau vs fo/f
figure('Name','Tau vs fo/f');
plot(FoF, Tau, 'k-', 'LineWidth', 3); hold on
yline(0,'k--','LineWidth',1);
grid on; grid minor
xlabel('$f_o/f$','Interpreter','latex','FontWeight','bold');
ylabel('$\tau=\frac{t_c-t_v}{T}$','Interpreter','latex','FontWeight','bold');
title('\tau vs $f_o/f$  (T = 1/f_{stylus})','Interpreter','latex');
set(gca,'LineWidth',1,'FontSize',16);

% Fig B — Lambda vs fo/f
figure('Name','Lambda vs fo/f');
plot(FoF, Lambda, 'Color',[0 0 0], 'LineWidth', 3); hold on
grid on; grid minor
xlabel('$f_o/f$','Interpreter','latex','FontWeight','bold');
ylabel('$\Lambda=V_{d,\max}/V_{s,\max}$','Interpreter','latex','FontWeight','bold');
title('\Lambda vs $f_o/f$','Interpreter','latex');
set(gca,'LineWidth',1,'FontSize',16);

% Fig C — Alpha (PRA) vs fo/f
figure('Name','Alpha (PRA) vs fo/f');
plot(FoF, Alpha_PRA, 'k-', 'LineWidth', 3); grid on; grid minor; hold on
xlabel('$f_o/f$','Interpreter','latex','FontWeight','bold');
ylabel('$\alpha=\left(V_e/V_p^\*\right)^2$','Interpreter','latex','FontWeight','bold');
title('Ejection efficiency $\alpha$ (PRA) vs $f_o/f$','Interpreter','latex');
set(gca,'LineWidth',1,'FontSize',16);

% (선택) CSV 저장 (기존)
tau_tbl = table(FoF(:), Tau(:), 'VariableNames', {'fo_over_f','tau'});
writetable(tau_tbl, 'tau_vs_foverf.csv');

lambda_tbl = table(FoF(:), Lambda(:), 'VariableNames', {'fo_over_f','Lambda'});
writetable(lambda_tbl, 'lambda_vs_foverf.csv');

% CSV 저장 (추가: PRA alpha)
alpha_pra_tbl = table(FoF(:), Alpha_PRA, ...
    'VariableNames', {'fo_over_f','alpha_pra'});
writetable(alpha_pra_tbl, 'alpha_pra_vs_foverf.csv');

% ---------- (추가) Normalize work and plot ----------
if ~isempty(Wraw)
    Wnorm = Wraw ./ max(abs(Wraw));   % sign 유지 정규화

    % Fig D — normalized energy transfer vs fo/f
    figure('Name','Normalized Work vs fo/f');
    plot(FoF, Wnorm, 'k-', 'LineWidth', 3); hold on
    yline(0,'k--','LineWidth',1);
    grid on; grid minor
    xlabel('$f_o/f$','Interpreter','latex','FontWeight','bold');
    ylabel('$\tilde{W}_{0\rightarrow t_e}$','Interpreter','latex','FontWeight','bold');
    title('Normalized energy transfer $\tilde{W}$ vs $f_o/f$','Interpreter','latex');
    set(gca,'LineWidth',1,'FontSize',16);

    % Fig E — Raw work vs fo/f
    figure('Name','Raw Work vs fo/f');
    plot(FoF, Wraw, 'k-', 'LineWidth', 3); hold on
    yline(0,'k--','LineWidth',1);
    grid on; grid minor
    xlabel('$f_o/f$','Interpreter','latex','FontWeight','bold');
    ylabel('$W_{0\rightarrow t_e}$','Interpreter','latex','FontWeight','bold');
    title('Raw work $W$ vs $f_o/f$','Interpreter','latex');
    set(gca,'LineWidth',1,'FontSize',16);

    % CSV 저장 (추가)
    work_tbl = table(FoF(:), Wraw(:), Wnorm(:), ...
        'VariableNames', {'fo_over_f','W_raw','W_norm'});
    writetable(work_tbl, 'work_vs_foverf.csv');
end

% =================== ODE function ===================
function ydot = supertwospring(t, y, pars)
% y = [x_s; v_s; x_d; v_d]
xs = y(1); vs = y(2);
xd = y(3); vd = y(4);

% Springs
Fk_s = pars.kstylus * (xs);             % stylus to ground
Fk_d = pars.kdroplet * (xd - xs);       % droplet to stylus

% Damping (relative)
Fd_d = pars.c * (vd - vs);

% EOM
axs = (-Fk_s - Fk_d - Fd_d) / pars.mstylus;
axd = ( -Fk_d - Fd_d )     / pars.mdroplet;

ydot = [vs; axs; vd; axd];
end
