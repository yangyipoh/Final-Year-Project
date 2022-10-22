clear all; close all; clc;

addpath('Data\')        % --------------- Create a folder for the data --------------
addpath('functions')

data = edfread('test3\test5.EDF'); % ------ Read .EDF file (in Data folder) ---------
fs = 1000;
ecg = cell2mat(data.ECG);

start_idx = 14089;      % ------------------ start_idx and end_idx ------------------
end_idx = 124579;
% start_idx = 1;
% end_idx = length(ecg);
ecg = ecg(start_idx:end_idx);

% wavelet filtering
[WT, F] = cwt(ecg, fs);
ecg_filt = icwt(WT, [], F, [0.5 150], 'SignalMean', mean(ecg));

% normalisation
range = max(ecg_filt(:)) - min(ecg_filt(:));
ecg_filt = (ecg_filt - min(ecg_filt(:))) / range;

% find R peaks
tm = 1:length(ecg);
tm = tm/fs;
[qrspeaks,locs] = findpeaks(ecg_filt,tm,'MinPeakHeight',0.85,...
            'MinPeakDistance',0.150);
qrspeaks = qrspeaks';

% ecg plot
figure(1)
plot(tm, ecg_filt)
hold on
plot(locs,qrspeaks,'ro')
xlim([0, 10])


% EDR --> RSA was found to be one of the best, but other methods can be tried
% EDR_physionet=edr1(0,ecg_filt,locs,fs);
% EDR_ramp=Rpeak_EDR(ecg_filt,locs,fs);
EDR_RSA=RSA_resp(ecg_filt,locs,fs);
% EDR_KPCA=KPCA_EDR(ecg_filt,locs,fs);

% EDR plot
figure(2)
% plot(EDR_physionet(:,1),EDR_physionet(:,2)),hold on,
% plot(EDR_ramp(:,1),EDR_ramp(:,2),'r'), hold on
plot(EDR_RSA(:,1),EDR_RSA(:,2),'r'), hold on
% plot(EDR_KPCA(:,1),EDR_KPCA(:,2),'k'),
legend('EDR-KPCA')
xlim([0 30])