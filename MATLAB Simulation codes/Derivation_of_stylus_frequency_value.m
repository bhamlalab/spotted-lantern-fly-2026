%% === Step-like theta(t) with t in MILLISECONDS → ERF fit → fs (Hz) ===
% t_ms: time in MILLISECONDS, theta: angle (deg or rad; frequency unaffected)
t_ms = [ ...
231.0555556
231.6111111
232.1666667
232.7222222
233.2777778
233.8333333
234.3888889
234.9444444
235.5
236.0555556
236.6111111
237.1666667
237.7222222
238.2777778
238.8333333
239.9444444
241.0555556
241.6111111
242.1666667
242.4444444
242.7222222
242.8333333
242.9444444
243.0555556
243.1666667
243.2222222
243.2777778
243.3333333
243.3888889
243.4444444
243.5
243.5555556
243.6111111
243.6666667
243.7222222
243.7777778
243.8333333
243.8888889
243.9444444
244
244.0555556
244.1111111
244.1666667
244.2222222
244.2777778
244.3333333
244.3888889
244.4444444
244.5
244.5555556
244.6111111
244.6666667
244.7222222
244.7777778
244.8333333
244.8888889
244.9444444
245
245.0555556
245.1111111
245.1666667
245.2222222
245.2777778
245.3333333
245.3888889
245.4444444
245.5
245.5555556
245.6111111
245.6666667
245.7222222
247.1666667
248.2777778
249.3888889
251.6111111
252.7222222
253.8333333
254.9444444];

theta = [ ...
85.656
87.31
87.353
88.094
87.077
87.385
88.235
88.938
89.919
90.4
91.407
91.532
92.631
92.932
93.007
92.903
92.922
90.671
89.131
88.029
84.374
82.968
81.137
78.688
76.838
74.683
72.358
70.903
68.235
66.001
63.701
59.187
55.781
52.778
47.751
43.569
39.361
35.573
30.975
25.312
22.109
19.175
16.773
12.163
10.278
7
5.973
4.358
3.121
1.785
1.267
0.402
-0.586
-2.22
-2.315
-4.81
-4.884
-6.305
-7.015
-9.509
-10.343
-10.867
-12.079
-12.213
-13.257
-12.951
-13.057
-13.152
-12.77
-12.834
-12.971
-12.649
-14.779
-14.915
-14.829
-15.367
-16.589
-16.641];

%% 0) Convert ms → s, sort, shift start time to 0
t_ms  = t_ms(:); theta = theta(:);
[ t_ms, ix ] = sort(t_ms); theta = theta(ix);
t = (t_ms - t_ms(1)) / 1000;         % *** milliseconds → seconds ***

%% 1) Detect transition window robustly (use slope magnitude)
g  = abs(gradient(theta, t));         % |dtheta/dt| with seconds
th = 0.2*max(g);                      % 20% of max slope
idx = find(g >= th);
if isempty(idx)
    error('Transition region not detected (slope too small).');
end
t10   = t(idx(1));
t90   = t(idx(end));
t0    = 0.5*(t10 + t90);
riseT = max(t90 - t10, eps);
mask  = (t >= t10 - 2*riseT) & (t <= t90 + 2*riseT);

%% 2) ERF fit: theta(t) = a*erf(b*t + c) + d
theta_hi = median(theta( max(1,idx(1)-5) : idx(1) ));
theta_lo = median(theta( idx(end) : min(numel(theta), idx(end)+5) ));
a0 = (theta_lo - theta_hi)/2;   % negative for downward step
d0 = (theta_hi + theta_lo)/2;
b0 = 4/max(riseT, eps);         % rough slope scale (1/s)
c0 = -b0*t0;

ft  = fittype('a*erf(b*x + c) + d', 'independent','x', ...
              'coefficients',{'a','b','c','d'});
opts = fitoptions('Method','NonlinearLeastSquares','Robust','LAR', ...
                  'StartPoint',[a0 b0 c0 d0]);

[mdl, gof] = fit(t(mask), theta(mask), ft, opts);
a = mdl.a; b = mdl.b; c = mdl.c; d = mdl.d;

%% 3) Analytic derivatives (smooth, noise-free)
td   = linspace(min(t(mask)), max(t(mask)), 2000).';
u    = b*td + c;
theta_fit = a*erf(u) + d;
omega_fit = (2*a*b/sqrt(pi)) * exp(-(u.^2));                  % dtheta/dt
alpha_fit = -(4*a*b^2/sqrt(pi)) * u .* exp(-(u.^2));          % d^2theta/dt^2

%% 4) Frequency from acceleration peak-to-peak (half-period)
[pkp, lp] = findpeaks(alpha_fit, td);              % +peaks (time positions lp)
[pkn, ln] = findpeaks(-alpha_fit, td); pkn = -pkn; % -peaks (time positions ln)
if isempty(pkp) || isempty(ln)
    error('Could not find sufficient peaks on acceleration.');
end
% pick the strongest +/- peaks
[~, ip] = max(pkp);   t_pos = lp(ip);
[~, in] = max(abs(pkn)); t_neg = ln(in);

Delta_t = abs(t_pos - t_neg);   % half period in SECONDS
fs_pp   = 1/(2*Delta_t);        % *** Hz ***

%% 5) (Optional) sanity-check: sine fit on acceleration
ft_s  = fittype('A*sin(2*pi*f*x + phi) + D', 'independent','x', ...
                'coefficients',{'A','f','phi','D'});
A0s = 0.5*(max(alpha_fit)-min(alpha_fit)); D0s = mean(alpha_fit);
opts_s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR', ...
                    'StartPoint',[A0s, fs_pp, 0, D0s]);
fs_sine = NaN; rmse_s = NaN;
try
    [mdl_s, gof_s] = fit(td, alpha_fit, ft_s, opts_s);
    fs_sine = mdl_s.f;
    if isfield(gof_s,'rmse'), rmse_s = gof_s.rmse; end
catch
    % keep NaN defaults
end

%% 6) Report (Hz)
fprintf('ERF fit RMSE                  = %.4g\n', gof.rmse);
fprintf('fs (alpha peak-to-peak)       = %.6f Hz\n', fs_pp);
fprintf('fs (sine fit on acceleration) = %.6f Hz (RMSE=%.3g)\n', fs_sine, rmse_s);

%% 7) Plots
figure('Name','ERF fit & derivatives (t in ms → converted to s)');
subplot(3,1,1);
plot(t, theta, 'k.', td, theta_fit, 'r-', 'LineWidth', 1.5); grid on
legend('\theta data','ERF fit','Location','best');
ylabel('\theta(t)'); title('Step (ERF) fit');

subplot(3,1,2);
plot(td, omega_fit, 'b-', 'LineWidth', 1.5); grid on
ylabel('\omega(t) = d\theta/dt');

subplot(3,1,3);
ypos = interp1(td, alpha_fit, t_pos);
yneg = interp1(td, alpha_fit, t_neg);
plot(td, alpha_fit, 'm-', 'LineWidth', 1.5); hold on
plot(t_pos, ypos, 'ro', t_neg, yneg, 'bo');
grid on; xlabel('time (s)'); ylabel('\alpha(t) = d^2\theta/dt^2');
title(sprintf('f_s \\approx %.4f Hz (from \\alpha peak-to-peak)', fs_pp));
