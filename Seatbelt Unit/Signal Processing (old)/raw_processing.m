clear variables, close all; clc;

%% read data
acc = readmatrix('test_breath_4.csv');
fs = 50;

%% Preprocessing
% smooting accelerometer
gx = 0;
gy = 0;
gz = 0;
alpha = 0.7;
linear_acc = zeros(size(acc));
for i = 1:size(linear_acc, 1)
    gx = alpha*gx + (1-alpha)*acc(i, 1);
    gy = alpha*gy + (1-alpha)*acc(i, 2);
    gz = alpha*gz + (1-alpha)*acc(i, 3);

    linear_acc(i, 1) = acc(i, 1) - gx;
    linear_acc(i, 2) = acc(i, 2) - gy;
    linear_acc(i, 3) = acc(i, 3) - gz;
end

mag_lin_acc = sqrt(linear_acc(:, 1).^2 + linear_acc(:, 2).^2 + linear_acc(:, 3).^2);
mag_lin_acc = mag_lin_acc(25:end);

t = 1:length(mag_lin_acc);
t = t/fs;
figure(1)
plot(t, mag_lin_acc)
title('Raw data')


order = 5;
wn = 2.5;
wn = wn/(fs/2);
[b, a] = butter(order, wn, 'low');
x_lp = filter(b, a, mag_lin_acc);

x_smooth = smooth(x_lp, 10);


figure(2)
plot(t, x_smooth)
title('Smooth data')

figure(3)
x_smooth = x_smooth - mean(x_smooth);
[S,F,T] = stft(x_smooth,fs,'Window',hamming(1024,'periodic'),'OverlapLength',700,'FrequencyRange','onesided');
waterfall(F,T,abs(S(:,:,1))')
xlim([0 3])
