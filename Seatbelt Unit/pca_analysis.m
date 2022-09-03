clear all, close all; clc;

%% read data
acc = readmatrix('test9.txt');
fs = 50;
acc = acc(:, 1:3);
%% Preprocessing
% Downsample to 10 Hz (factor = 5)
% lp filter
order = 8;
wn = 2.0;
wn = wn/(fs/2);
[b, a] = butter(order, wn, 'low');
acc_lp = filter(b, a, acc, [], 1);

% downsample
factor = 5;
acc_ds = downsample(acc_lp, factor);
fs = fs/factor;

% smooting accelerometer
gx = 0;
gy = 0;
gz = 0;
alpha = 0.8;
tm = 0.06;
linear_acc = zeros(size(acc_ds));
for i = 1:size(linear_acc, 1)
    gx = alpha*gx + (1-alpha)*acc_ds(i, 1);
    gy = alpha*gy + (1-alpha)*acc_ds(i, 2);
    gz = alpha*gz + (1-alpha)*acc_ds(i, 3);

    linear_acc(i, 1) = acc_ds(i, 1) - gx;
    linear_acc(i, 2) = acc_ds(i, 2) - gy;
    linear_acc(i, 3) = acc_ds(i, 3) - gz;
end
% linear_acc = linear_acc(round(delay*fs):end, :);     % getting rid of the first few samples
linear_acc(abs(linear_acc) > tm) = NaN; 


coeff = pca(linear_acc);
breathing = coeff(:, 1)'*linear_acc';
breathing = fillmissing(breathing, "linear");

order = 5;
wn = 0.8;
wn = wn/(fs/2);
[b, a] = butter(order, wn, 'low');
breathing = filter(b, a, breathing);

delay = 3.0;
breathing = breathing(round(delay*fs):end);

t = 1:size(acc, 1);
t = t/50;
figure(1)
sgtitle('Raw data')
subplot(3, 1, 1)
plot(t, acc(:, 1))
title('x axis')
subplot(3, 1, 2)
plot(t, acc(:, 2))
title('y axis')
subplot(3, 1, 3)
plot(t, acc(:, 3))
title('z axis')

t = 1:length(breathing);
t = t/fs;
figure(2)
plot(t, breathing)
title('Dominant axis')

figure(3)
[S,F,T] = stft(breathing,fs,'Window',hamming(64,'periodic'),'OverlapLength',60,'FrequencyRange','onesided');
waterfall(F,T,abs(S(:,:,1))')
xlim([0 3])
title('STFT')

%% Calculate breathing rate
S_mag = abs(S);

[~, idx] = max(S_mag, [], 1);
for i = 1:length(idx)
    fprintf('Time: %.2fs, Breathing rate: %.2fHz\n', T(i), F(idx(i)))
end

% Seat belt --> Adaptive filtering
