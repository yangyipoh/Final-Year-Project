clear; close all; clc;

trial = 3;
test = 4;  % 3,4=moving, 2,5=stationary, bumps=3,5

% start_idx = 3132;  % third
% end_idx = 59198;
% start_idx = 4753;  % second
% end_idx = 51036;
% start_idx = 13941;
% end_idx = 83699;   % first
% start_idx = 4787;
% end_idx = 30296;   % fifth
% start_idx = 9511;
% end_idx = 53829;    % six
% start_idx = 26344;
% end_idx = 67287;   % seven
% start_idx = 9671;
% end_idx = 45727;    % eight
% start_idx = 57416;
% end_idx = 108080;    % nine
% start_idx = 15646;
% end_idx = 49280;    % one
% start_idx = 13794;
% end_idx = 60931;    % one
% start_idx = 9697;
% end_idx = 56987;   % two
% start_idx = 5104;
% end_idx = 55905;   % three
start_idx = 8825;
end_idx = 54677;   % four
% start_idx = 15638;
% end_idx = 123769;   %five

%% Add external path
addpath(sprintf("Data\\test%d\\", trial));     % select trial
addpath("functions\");      % adaptive filtering functions

%% ACCELEROMETER -- Read data
fname = sprintf("test%d.txt", test);
data = readmatrix(fname);
fs_resp = 50;

order = 4;
wn = 1.0;
wn = wn/(fs_resp/2);
[b, a] = butter(order, wn, 'low');
data = filter(b, a, data, [], 1);
% data = data(50*8:end, :);

% downsample
factor = 5;
data = downsample(data, factor);
fs_resp = fs_resp/factor;

data_main = data(:, 1:3);
data_ref = data(:, 4:6);

%% Processing

% analyse first 20 seconds
% offset = 25*50;
% data_main = data_main((1+offset):(50*30+offset), :);
% data_ref = data_ref((1+offset):(50*30+offset), :);

% align axis
[n, ~] = size(data_main);
[rot_matrix, t] = calib_rotation(data_main', data_ref');
data_main = (rot_matrix*data_main') + repmat(t, 1, n);
data_main = data_main';

% smooth data
% data_main = smoothdata(data_main, 1, 'gaussian', 10);
% data_ref = smoothdata(data_ref, 1, 'gaussian', 10);

%% adaptive filtering
[n, ~] = size(data_main);
filtered_acc = zeros(n, 3);

for i = 1:3
%     e = data_main(:, i) - data_ref(:, i);

    weight_L = 5;
    RLS = dsp.RLSFilter('Length', weight_L, 'Method','Conventional RLS');
    [y, e] = RLS(data_ref(:, i),data_main(:, i));

%     mu = 0.07;                    % Step size
%     po = 8;                      % Projection order
%     offset = 0.01;               % Offset for covariance matrix
%     apf = dsp.AffineProjectionFilter('Length', 32, ...
%         'StepSize', mu, 'ProjectionOrder', po, ...
%         'InitialOffsetCovariance',offset);
%     [y,e] = apf(data_ref(:, i), data_main(:, i));
    
%     lam = 0.995;
%     del = 1.2;
%     alf = dsp.AdaptiveLatticeFilter('Length', 32, ...
%     'ForgettingFactor', lam, 'InitialPredictionErrorPower', del);
%     [y,e] = alf(data_ref(:, i), data_main(:, i));
    
%     mu = 0.008;
%     wn = 2.0;
%     wn = wn/(fs_resp/2);
%     b  = fir1(5,wn);
%     fxlms = dsp.FilteredXLMSFilter(32, 'StepSize', mu, 'LeakageFactor', ...
%      1, 'SecondaryPathCoefficients', b);
%     [y,e] = fxlms(data_ref(:, i), data_main(:, i));

    filtered_acc(:, i) = e;
%     filtered_acc(:, i) = data_main(:, i);
end

gx = 0;
gy = 0;
gz = 0;
alpha = 0.85;
linear_acc = zeros(size(filtered_acc));
for i = 1:size(linear_acc, 1)
    gx = alpha*gx + (1-alpha)*filtered_acc(i, 1);
    gy = alpha*gy + (1-alpha)*filtered_acc(i, 2);
    gz = alpha*gz + (1-alpha)*filtered_acc(i, 3);

    linear_acc(i, 1) = filtered_acc(i, 1) - gx;
    linear_acc(i, 2) = filtered_acc(i, 2) - gy;
    linear_acc(i, 3) = filtered_acc(i, 3) - gz;
end
order = 5;
wn = 0.7;
wn = wn/(fs_resp/2);
[b, a] = butter(order, wn, 'low');
linear_acc = filter(b, a, linear_acc);
linear_acc = linear_acc(20:end, :);

linear_acc_plot = linear_acc;
tm_movement = 0.06;
linear_acc(abs(linear_acc) > tm_movement) = NaN;

coeff = pca(linear_acc);
breathing = coeff(:, 1)'*linear_acc';
breathing = fillmissing(breathing, "linear");

% order = 5;
% wn = 0.28;
% wn = wn/(fs_resp/2);
% [b, a] = butter(order, wn, 'low');
% breathing = filter(b, a, breathing);

% delay = 1.0;
% breathing = breathing(round(delay*fs_resp):end);

%% Wavelet decomposition

% decompose signal using discrete wavelet transform
[WT, F] = cwt(breathing, fs_resp);
breathing = icwt(WT, [], F, [0.15 0.43], 'SignalMean', mean(breathing));

breathing = [zeros(1, 19)*breathing(1) breathing];
%% EDR
fname = sprintf("test%d.EDF", test);
data = edfread(fname);
fs_edr = 1000;
ecg = cell2mat(data.ECG);
ecg = ecg(start_idx:end_idx);

[WT, F] = cwt(ecg, fs_edr);
ecg = icwt(WT, [], F, [0.5 150], 'SignalMean', mean(ecg));

tm = 1:length(ecg);
tm = tm/fs_edr;

range = max(ecg(:)) - min(ecg(:));
ecg = (ecg - min(ecg(:))) / range;

% locate R peaks
[qrspeaks,locs_breathing] = findpeaks(ecg,tm,'MinPeakHeight',0.85,...
            'MinPeakDistance',0.150);
qrspeaks = qrspeaks';

% EDR_ramp=Rpeak_EDR(ecg,locs,fs_edr);
% EDR_KPCA=KPCA_EDR(ecg,locs,fs_edr);
EDR_RSA=RSA_resp(ecg,locs_breathing,fs_edr);

%% Plot

figure(1)
plot(tm, ecg)
hold on
plot(locs_breathing,qrspeaks,'ro')

t = 1:size(data_main, 1);
t = t/10;
figure(2)
subplot(3, 1, 1)
plot(t, data_main(:, 1))
hold on
plot(t, data_ref(:, 1))
xlabel('time (s)')
ylabel('acceleration (g)')
title('x-axis')
% xlim([35 45])

subplot(3, 1, 2)
plot(t, data_main(:, 2))
hold on
plot(t, data_ref(:, 2))
xlabel('time (s)')
ylabel('acceleration (g)')
title('y-axis')
% xlim([35 45])

subplot(3, 1, 3)
plot(t, data_main(:, 3))
hold on
plot(t, data_ref(:, 3))
xlabel('time (s)')
ylabel('acceleration (g)')
title('z-axis')
% xlim([35 45])
sgtitle('Raw Accelerometer data')


t = 1:size(linear_acc, 1);
t = t/10;
figure(3)
sgtitle('Noise cancelled data')
subplot(3, 1, 1)
plot(t, linear_acc(:, 1))
title('x axis')
xlabel('time (s)')
ylabel('acceleration (g)')
% xlim([35 45])
subplot(3, 1, 2)
plot(t, linear_acc(:, 2))
title('y axis')
xlabel('time (s)')
ylabel('acceleration (g)')
% xlim([35 45])
subplot(3, 1, 3)
plot(t, linear_acc(:, 3))
title('z axis')
xlabel('time (s)')
ylabel('acceleration (g)')
% xlim([35 45])

t = 1:length(breathing);
t = t/fs_resp;
figure(4)
yyaxis left
plot(t, -breathing)
ylabel('Acceleration')

yyaxis right
plot(EDR_RSA(:,1),EDR_RSA(:,2),'r'), hold on
% plot(EDR_KPCA(:,1),2*EDR_KPCA(:,2),'k')
ylabel('voltage')
title('Respiration waveforms')
% % xlim([0 30])
xlabel('time')
% xlim([35 45])
%% Calculate breathing rate
[pks, locs_breathing] = findpeaks(-breathing, 'threshold', 0);
yyaxis left
hold on
plot(t(locs_breathing), pks, '*')
locs_breathing = locs_breathing/fs_resp;

[pks, locs_actual] = findpeaks(EDR_RSA(:, 2), 'threshold', 0);
yyaxis right
hold on
plot(t(locs_actual), pks, '*')
fs_edr = 1/mean(diff(EDR_RSA(:, 1)));
locs_actual = locs_actual/fs_edr;
legend('Accelerometer', 'Accelerometer peaks','EDR-RSA', 'EDR peaks')

total_time = t(end);
plt_res = zeros(2, 0);
idx = 1;
window_size = 12;
step_size = 1;
for i = window_size:step_size:total_time
    exp_res = locs_breathing >= i-window_size & locs_breathing <= i;
    exp_res = mean(diff(locs_breathing(exp_res)));
    act_res = locs_actual >= i-window_size & locs_actual <= i;
    act_res = mean(diff(locs_actual(act_res)));
    plt_res(:, idx) = [1/exp_res*60; 1/act_res*60];
    idx = idx + 1;
end

figure(5)
plot(plt_res(1, :), plt_res(2, :), '*', 'MarkerSize', 10)
hold on
plot([15 30], [15 30])
xlabel('Accelerometer breathing rate')
ylabel('EDR breathing rate')
title('Breathing rate correlation')


average_val = (plt_res(1, :) + plt_res(2, :))*0.5;
diff_val = plt_res(1, :) - plt_res(2, :);
avg_diff = mean(diff_val);
std_lower = avg_diff - 1.96*std(diff_val);
std_upper = avg_diff + 1.96*std(diff_val);
figure(6)
scatter(average_val, diff_val)
hold on
yline(avg_diff, 'Color', 'k')
yline(std_lower, 'Color', 'r')
yline(std_upper, 'Color', 'r')
xlabel('Average')
ylabel('Difference')
title('Bland-Altman plot of Breathing rate')