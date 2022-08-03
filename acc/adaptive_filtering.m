clear all; close all; clc;

addpath("data\test1\");
addpath("functions\");
%% Align axis
% calib_data = readmatrix('calib.txt');
% calib_main = calib_data(:, 1:3)';
% calib_ref = calib_data(:, 4:6)';

%% Read data
data = readmatrix('test4.txt');
data_main = data(:, 1:3);
data_ref = data(:, 4:6);
fs = 50;

%% Processing
% align axis
[n, ~] = size(data_main);
[rot_matrix, t] = calib_rotation(data_main', data_ref');
data_main = (rot_matrix*data_main') + repmat(t, 1, n);
data_main = data_main';
% data_ref = data_ref';

% smooth data
data_main = smoothdata(data_main, 1, 'gaussian', 25);
data_ref = smoothdata(data_ref, 1, 'gaussian', 25);

figure(4)
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

% RLS adaptive filtering
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
%     del = 1;
%     alf = dsp.AdaptiveLatticeFilter('Length', 32, ...
%     'ForgettingFactor', lam, 'InitialPredictionErrorPower', del);
%     [y,e] = alf(data_ref(:, i), data_main(:, i));
    
    mu = 0.008;
    wn = 2.0;
    wn = wn/(fs/2);
    b  = fir1(5,wn);
    fxlms = dsp.FilteredXLMSFilter(32, 'StepSize', mu, 'LeakageFactor', ...
     1, 'SecondaryPathCoefficients', b);
    [y,e] = fxlms(data_ref(:, i), data_main(:, i));

    filtered_acc(:, i) = e;
end

% downsample to 10 Hz
order = 10;
wn = 2.0;
wn = wn/(fs/2);
[b, a] = butter(order, wn, 'low');
filtered_acc_lp = filter(b, a, filtered_acc, [], 1);
factor = 5;
filtered_acc_ds = downsample(filtered_acc_lp, factor);
fs = fs/factor;


tm = 0.05;
filtered_acc_ds(abs(filtered_acc_ds) > tm) = NaN;

linear_acc = filtered_acc_ds;

coeff = pca(linear_acc);
breathing = coeff(:, 1)'*linear_acc';
breathing = fillmissing(breathing, "spline");

order = 5;
wn = 0.5;
wn = wn/(fs/2);
[b, a] = butter(order, wn, 'low');
breathing = filter(b, a, breathing);

delay = 5.0;
breathing = breathing(round(delay*fs):end);

t = 1:size(filtered_acc, 1);
t = t/10;
figure(1)
sgtitle('Raw data')
subplot(3, 1, 1)
plot(t, filtered_acc(:, 1))
title('x axis')
ylim([-tm tm])
subplot(3, 1, 2)
plot(t, filtered_acc(:, 2))
title('y axis')
ylim([-tm tm])
subplot(3, 1, 3)
plot(t, filtered_acc(:, 3))
title('z axis')
ylim([-tm tm])

t = 1:length(breathing);
t = t/fs;
figure(2)
plot(t, breathing)
title('Dominant axis')
xlim([0 30])

figure(3)
[S,F,T] = stft(breathing,fs,'Window',hamming(128,'periodic'),'OverlapLength',100,'FrequencyRange','onesided');
waterfall(F,T,abs(S(:,:,1))')
xlim([0 3])
title('STFT')

%% Calculate breathing rate
S_mag = abs(S);
S_mag(1, :) = 0;
[~, idx] = max(S_mag, [], 1);
for i = 1:length(idx)
    fprintf('Time: %.2fs, Breathing rate: %.2f bpm\n', T(i), 60*F(idx(i)))
end

