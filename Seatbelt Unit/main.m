clear; close all; clc;

trial = 2;
test = 1;

%% Add external path
addpath(sprintf("Data\\test%d\\", trial));     % select trial
addpath("functions\");      % adaptive filtering functions
addpath('edrfunctions\');   % ECG derived respiration functions

%% ACCELEROMETER -- Read data
fname = sprintf("test%d.txt", test);
data = readmatrix(fname);
data_main = data(:, 1:3);
data_ref = data(:, 4:6);
fs_resp = 50;

%% Processing
% truncate the first and last data
data_main = data_main(100:end-100, :);
data_ref = data_ref(100:end-100, :);

% analyse first 20 seconds
offset = 25*50;
% data_main = data_main((1+offset):(50*30+offset), :);
% data_ref = data_ref((1+offset):(50*30+offset), :);

% align axis
[n, ~] = size(data_main);
[rot_matrix, t] = calib_rotation(data_main', data_ref');
data_main = (rot_matrix*data_main') + repmat(t, 1, n);
data_main = data_main';

% smooth data
data_main = smoothdata(data_main, 1, 'gaussian', 25);
data_ref = smoothdata(data_ref, 1, 'gaussian', 25);

%% adaptive filtering
[n, ~] = size(data_main);
filtered_acc = zeros(n, 3);

for i = 1:3
%     e = data_main(:, i) - data_ref(:, i);

%     weight_L = 10;
%     RLS = dsp.RLSFilter('Length', weight_L, 'Method','Conventional RLS');
%     [y, e] = RLS(data_ref(:, i),data_main(:, i));

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
    
    mu = 0.008;
    wn = 2.0;
    wn = wn/(fs_resp/2);
    b  = fir1(5,wn);
    fxlms = dsp.FilteredXLMSFilter(32, 'StepSize', mu, 'LeakageFactor', ...
     1, 'SecondaryPathCoefficients', b);
    [y,e] = fxlms(data_ref(:, i), data_main(:, i));

    filtered_acc(:, i) = e;
end

% downsample to 10 Hz
order = 10;
wn = 2.0;
wn = wn/(fs_resp/2);
[b, a] = butter(order, wn, 'low');
filtered_acc_lp = filter(b, a, filtered_acc, [], 1);
factor = 5;
filtered_acc_ds = downsample(filtered_acc_lp, factor);
fs_resp = fs_resp/factor;


tm_movement = 0.03;
filtered_acc_ds(abs(filtered_acc_ds) > tm_movement) = NaN;

linear_acc = filtered_acc_ds;

coeff = pca(linear_acc);
breathing = coeff(:, 1)'*linear_acc';
breathing = fillmissing(breathing, "linear");

order = 5;
wn = 0.28;
wn = wn/(fs_resp/2);
[b, a] = butter(order, wn, 'low');
breathing = filter(b, a, breathing);

% delay = 1.0;
% breathing = breathing(round(delay*fs_resp):end);

%% Wavelet decomposition

% decompose signal using discrete wavelet transform
[WT, F] = cwt(breathing, fs_resp);
breathing = icwt(WT, [], F, [0.16 0.30], 'SignalMean', mean(breathing));


%% EDR
fname = sprintf("test%d.EDF", test);
data = edfread(fname);
fs_edr = 1000;
ecg = cell2mat(data.ECG);
% start_idx = 3132;  % third
% end_idx = 59198;
% start_idx = 4753;  % second
% end_idx = 51036;
start_idx = 13941;
end_idx = 83699;   % first
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
EDR_KPCA=KPCA_EDR(ecg,locs,fs_edr);
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

t = 1:size(filtered_acc, 1);
t = t/50;
figure(3)
sgtitle('Raw data')
subplot(3, 1, 1)
plot(t, filtered_acc(:, 1))
title('x axis')
% ylim([-tm_movement tm_movement])
subplot(3, 1, 2)
plot(t, filtered_acc(:, 2))
title('y axis')
% ylim([-tm_movement tm_movement])
subplot(3, 1, 3)
plot(t, filtered_acc(:, 3))
title('z axis')
% ylim([-tm_movement tm_movement])
% 
t = 1:length(breathing);
t = t/fs_resp;
figure(4)
subplot(2, 1, 1)
plot(t, breathing)
title('Accelerometer')
xlim([0 30])

subplot(2, 1, 2)
plot(EDR_RSA(:,1),EDR_RSA(:,2),'r'), hold on
plot(EDR_KPCA(:,1),2*EDR_KPCA(:,2),'k')
legend('EDR-RSA', 'EDR-KPCA')
title('EDR')
xlim([0 30])

figure(5)
[S,F,T] = stft(breathing,fs_resp,'Window',hamming(256,'periodic'),'OverlapLength',240,'FrequencyRange','onesided');
waterfall(F,T,abs(S(:,:,1))')
xlim([0 1])
title('STFT')

%% Calculate breathing rate
S_mag = abs(S);
S_mag(1, :) = 0;
[~, idx] = max(S_mag, [], 1);
for i = 1:length(idx)-1
    fprintf('Time: %.2fs, Breathing rate: %.2f bpm\n', T(i), 60*F(idx(i)))
end

