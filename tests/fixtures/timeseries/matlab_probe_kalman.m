% matlab_probe_kalman.m — MATLAB R2018b ssm 批处理冒烟探针 (2026-08-20)
% 目的: ① 确认 -r 批处理可编程访问 ② ssm/estimate/filter/smooth 链路在场
%       ③ 留下确定性数值锚雏形 (rng 固定, %.17g 全精度) — v1.8 Kalman 对照生态预研
% 模型: 局部水平 y_t = a_t + e_t (D=Var(e)), a_t = a_{t-1} + u_t (B=Var(u))
diary('matlab_probe_kalman.log');
rng(42, 'twister');
T = 100;
y = cumsum(randn(T,1)) + 0.5*randn(T,1);
Mdl = ssm(1, NaN, 1, NaN);   % A=1, B=est, C=1, D=est
[estMdl, estParams, logL] = estimate(Mdl, y, [0.5; 1.0], 'Display', 'off');
fprintf('MATLAB_VERSION = %s\n', version);
fprintf('T = %d\n', T);
fprintf('EST_PARAMS = %.17g %.17g\n', estParams(1), estParams(2));
fprintf('LOGLIK = %.17g\n', logL);
[yF, ~] = filter(estMdl, y);      % Kalman filter (一步预测/滤波态)
fprintf('FILT_END = %.17g\n', yF(end));
[yS, ~] = smooth(estMdl, y);      % Kalman smoother (固定区间平滑)
fprintf('SMOOTH_END = %.17g\n', yS(end));
fprintf('Y_FIRST5 = %.17g %.17g %.17g %.17g %.17g\n', y(1:5));
fprintf('PROBE_OK\n');
diary off;
exit;
