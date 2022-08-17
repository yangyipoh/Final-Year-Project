clear all; close all; clc;

addpath('data\')
addpath('otherfunctions')

data = edfread('test6.EDF');
fs = 1000;
ecg = cell2mat(data.ECG);
ecg = ecg(19.11*fs:61*fs);

tm = 1:length(ecg);
tm = tm/fs;

range = max(ecg(:)) - min(ecg(:));
ecg = (ecg - min(ecg(:))) / range;

% locate R peaks
[qrspeaks,locs] = findpeaks(ecg,tm,'MinPeakHeight',0.8,...
            'MinPeakDistance',0.150);
qrspeaks = qrspeaks';

figure(1)
plot(tm, ecg)
hold on
plot(locs,qrspeaks,'ro')
% xlim([0, 2])


%% EDR
% EDR_physionet=edr1(0,ecg,locs,fs);
EDR_ramp=Rpeak_EDR(ecg,locs,fs);

figure(2)
% plot(EDR_physionet(:,1),EDR_physionet(:,2)),hold on, 
plot(EDR_ramp(:,1),EDR_ramp(:,2),'k'),
legend('physionet EDR','EDR-RAMP')
xlim([0 30])
