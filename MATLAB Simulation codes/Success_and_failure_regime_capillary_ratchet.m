%% =========================================================
%  Δκ (1/mm) vs θ (deg) with theory boundaries
%  + Export generated results to Excel
%
%  - Plots:
%     (1) Theory boundary curves (F_net = 0) for nymph and adult
%     (2) Δκ_min thresholds for nymph and adult
%     (3) Data points (nymph, adult)
%
%  - Exports to Excel:
%     Sheet 1: TheoryCurves   (dk_mm, thetaCrit_nymph_deg, thetaCrit_adult_deg)
%     Sheet 2: Thresholds     (dkmin_nymph_mm, dkmin_adult_mm)
%     Sheet 3: DataPoints     (group, dk_mm, theta_deg)
%% =========================================================

clear; clc; close all;

%% ---------- Physical/Model parameters ----------
gamma = 0.04;
mu    = 1.5e-3;
W     = 200e-6;
v     = 0.4;

%% ---------- Nymph parameters ----------
D_n   = 0.6e-3;
A_n   = pi/4 * D_n^2;
ld_n  = 0.2e-3;
x_n   = 0.5e-3;
CAH_n = 0.212;

% Nymph data point (from your note)
dk_obs_n_mm = 1.0;     % Δκ in 1/mm
theta_obs_n = 50;      % deg

%% ---------- Adult parameters ----------
D_a   = 0.7e-3;
A_a   = pi/4 * D_a^2;
ld_a  = 1.0e-3;
x_a   = 0.5e-3;
CAH_a = 0.27;

% Adult data point (from your note)
dk_obs_a_mm = 0.0;     % Δκ in 1/mm
theta_obs_a = 100;     % deg

%% ---------- Plot settings ----------
xMax_plot = 1.2;   % Δκ max (1/mm)
yMax_plot = 120;   % θ max (deg)
Ncurve    = 800;   % resolution of theory curves
epsFactor = 1.02;  % start slightly above threshold to avoid denom=0

%% ---------- Model functions ----------
dk_threshold_m = @(A, CAH) (W*CAH) / (2*A); % returns 1/m

theta_crit_rad = @(dk_m, A, ld, x, CAH) ...
    (mu*W*ld*v) ./ ( x .* (2*gamma*dk_m.*A - gamma*W*CAH) ); % rad

%% ---------- Thresholds (convert to 1/mm) ----------
dkmin_n_mm = dk_threshold_m(A_n, CAH_n) / 1e3; % 1/mm
dkmin_a_mm = dk_threshold_m(A_a, CAH_a) / 1e3; % 1/mm

%% ---------- Build theory curves ----------
dk_mm = linspace(0, xMax_plot, Ncurve).'; % column vector (1/mm)

% Preallocate with NaNs (below threshold is undefined, keep NaN)
thetaCrit_n_deg = nan(size(dk_mm));
thetaCrit_a_deg = nan(size(dk_mm));

% Nymph curve region (dk > dkmin)
mask_n = dk_mm >= dkmin_n_mm * epsFactor;
dk_n_m = dk_mm(mask_n) * 1e3; % 1/m
thetaCrit_n_deg(mask_n) = theta_crit_rad(dk_n_m, A_n, ld_n, x_n, CAH_n) * (180/pi);

% Adult curve region (dk > dkmin)
mask_a = dk_mm >= dkmin_a_mm * epsFactor;
dk_a_m = dk_mm(mask_a) * 1e3; % 1/m
thetaCrit_a_deg(mask_a) = theta_crit_rad(dk_a_m, A_a, ld_a, x_a, CAH_a) * (180/pi);

%% ---------- Plot ----------
figure('Position',[100 100 780 560]); hold on; box on;

plot(dk_mm, thetaCrit_n_deg, 'LineWidth', 1.6);
plot(dk_mm, thetaCrit_a_deg, 'LineWidth', 1.6);

scatter(dk_obs_n_mm, theta_obs_n, 70, 'filled');                 % nymph
scatter(dk_obs_a_mm, theta_obs_a, 80, 'filled', 'Marker','s');   % adult

xline(dkmin_n_mm, '--', sprintf('Nymph \\Delta\\kappa_{min}=%.3f 1/mm', dkmin_n_mm), 'LineWidth', 1);
xline(dkmin_a_mm, '--', sprintf('Adult \\Delta\\kappa_{min}=%.3f 1/mm', dkmin_a_mm), 'LineWidth', 1);

xlabel('Curvature difference \Delta\kappa (1/mm)');
ylabel('Stylus angle \theta (deg)');
xlim([0 xMax_plot]);
ylim([0 yMax_plot]);

legend({'Theory (nymph params)', 'Theory (adult params)', ...
        'Data: nymph', 'Data: adult', ...
        'Nymph \Delta\kappa_{min}', 'Adult \Delta\kappa_{min}'}, ...
        'Location','northeast');

title('\Delta\kappa vs \theta (theory + data)');
set(gca,'FontSize',11);

%% ---------- Export generated results to Excel ----------
% 1) Theory curves table
Tcurve = table(dk_mm, thetaCrit_n_deg, thetaCrit_a_deg, ...
    'VariableNames', {'dk_mm', 'thetaCrit_nymph_deg', 'thetaCrit_adult_deg'});

% 2) Thresholds table
Tth = table(dkmin_n_mm, dkmin_a_mm, ...
    'VariableNames', {'dkmin_nymph_mm', 'dkmin_adult_mm'});

% 3) Data points table
group = ["nymph"; "adult"];
dk_data_mm = [dk_obs_n_mm; dk_obs_a_mm];
theta_data_deg = [theta_obs_n; theta_obs_a];

Tdata = table(group, dk_data_mm, theta_data_deg, ...
    'VariableNames', {'group', 'dk_mm', 'theta_deg'});

% Choose save location
[saveFile, savePath] = uiputfile('ratchet_theory_export.xlsx', 'Save Excel file as');
if isequal(saveFile,0)
    disp("Export canceled.");
else
    outXlsx = fullfile(savePath, saveFile);

    writetable(Tcurve, outXlsx, 'Sheet', 'TheoryCurves');
    writetable(Tth,    outXlsx, 'Sheet', 'Thresholds');
    writetable(Tdata,  outXlsx, 'Sheet', 'DataPoints');

    fprintf("Exported Excel:\n  %s\n", outXlsx);
end

%% ---------- Quick console prints ----------
fprintf("Nymph Δκ_min = %.4f (1/mm)\n", dkmin_n_mm);
fprintf("Adult Δκ_min = %.4f (1/mm)\n", dkmin_a_mm);
