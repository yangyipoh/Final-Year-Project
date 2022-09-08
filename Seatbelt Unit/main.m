clear; close all; clc;

trial = 3;
test = 4;

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
% start_idx = 14136;
% end_idx = 124579;   %five

%% Add external path
addpath(sprintf("Data\\test%d\\", trial));     % select trial
addpath("functions\");      % adaptive filtering functions
addpath('edrfunctions\');   % ECG derived respiration functions

%% ACCELEROMETER -- Read data
fname = sprintf("test%d.txt", test);
data = readmatrix(fname);
fs_resp = 50;

order = 4;
wn = 1.0;
wn = wn/(fs_resp/2);
[b, a] = butter(order, wn, 'low');
data = filter(b, a, data, [], 1);
data = data(30:end, :);

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

    weight_L = 10;
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

breathing = [ones(1, 19)*breathing(1) breathing];
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
[qrspeaks,locs] = findpeaks(ecg,tm,'MinPeakHeight',0.85,...
            'MinPeakDistance',0.150);
qrspeaks = qrspeaks';

% EDR_ramp=Rpeak_EDR(ecg,locs,fs_edr);
% EDR_KPCA=KPCA_EDR(ecg,locs,fs_edr);
EDR_RSA=RSA_resp(ecg,locs,fs_edr);

%% Plot

figure(1)
plot(tm, ecg)
hold on
plot(locs,qrspeaks,'ro')

figure(2)
subplot(3, 1, 1)
plot(data_main(:, 1))
hold on
plot(data_ref(:, 1))

subplot(3, 1, 2)
plot(data_main(:, 2))
hold on
plot(data_ref(:, 2))

subplot(3, 1, 3)
plot(data_main(:, 3))
hold on
plot(data_ref(:, 3))

t = 1:size(linear_acc, 1);
t = t/10;
figure(3)
sgtitle('Raw data')
subplot(3, 1, 1)
plot(t, linear_acc(:, 1))
title('x axis')
% ylim([-tm_movement tm_movement])
subplot(3, 1, 2)
plot(t, linear_acc(:, 2))
title('y axis')
% ylim([-tm_movement tm_movement])
subplot(3, 1, 3)
plot(t, linear_acc(:, 3))
title('z axis')
% ylim([-tm_movement tm_movement])
% 
t = 1:length(breathing);
t = t/fs_resp;
figure(4)
yyaxis left
plot(t, -breathing)
title('Accelerometer')
% xlim([0 30])

yyaxis right
plot(EDR_RSA(:,1),EDR_RSA(:,2),'r'), hold on
% plot(EDR_KPCA(:,1),2*EDR_KPCA(:,2),'k')
legend('Accelerometer', 'EDR-RSA', 'EDR-KPCA')
title('EDR')
xlim([0 30])

% figure(5)
% [S,F,T] = stft(breathing,fs_resp,'Window',hamming(256,'periodic'),'OverlapLength',240,'FrequencyRange','onesided');
% waterfall(F,T,abs(S(:,:,1))')
% xlim([0 1])
% title('STFT')

%% Calculate breathing rate
% S_mag = abs(S);
% S_mag(1, :) = 0;
% [~, idx] = max(S_mag, [], 1);
% for i = 1:length(idx)
%     fprintf('Time: %.2fs, Breathing rate: %.2f bpm\n', T(i), 60*F(idx(i)))
% end

