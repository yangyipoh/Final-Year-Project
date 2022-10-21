clear all; close all; clc;

addpath('Data\test3\')
addpath('edrfunctions')

data = edfread('test5.EDF');
fs = 1000;
ecg = cell2mat(data.ECG);

%%
% start_idx = 14136;
% end_idx = 124579;
start_idx = 1;
end_idx = length(ecg);
ecg = ecg(start_idx:end_idx);

[WT, F] = cwt(ecg, fs);
ecg_filt = icwt(WT, [], F, [0.5 150], 'SignalMean', mean(ecg));
ecg_filt = ecg;

% ecg = detrend(ecg, 6);

tm = 1:length(ecg);
tm = tm/fs;

range = max(ecg_filt(:)) - min(ecg_filt(:));
ecg_filt = (ecg_filt - min(ecg_filt(:))) / range;

% locate R peaks
[qrspeaks,locs] = findpeaks(ecg_filt,tm,'MinPeakHeight',0.85,...
            'MinPeakDistance',0.150);
qrspeaks = qrspeaks';

figure(1)
plot(tm, ecg_filt)
hold on
plot(locs,qrspeaks,'ro')
xlim([0, 10])


%% EDR
% EDR_physionet=edr1(0,ecg_filt,locs,fs);
% EDR_ramp=Rpeak_EDR(ecg_filt,locs,fs);
EDR_RSA=RSA_resp(ecg_filt,locs,fs);
% EDR_KPCA=KPCA_EDR(ecg_filt,locs,fs);

figure(2)
% plot(EDR_physionet(:,1),EDR_physionet(:,2)),hold on,
% plot(EDR_ramp(:,1),EDR_ramp(:,2),'r'), hold on
plot(EDR_RSA(:,1),EDR_RSA(:,2),'y'), hold on
% plot(EDR_KPCA(:,1),EDR_KPCA(:,2),'k'),
legend('EDR-RSA', 'EDR-KPCA')
xlim([0 30])