load('test_data.mat')% Data contains ECG (ECG), r-peak locations (R) and RR intervals (RR),
addpath('otherfunctions')
% EDR code slightly modified from Physionet:
EDR_physionet=edr1(0,ECG,R,fs);
% R-peak ramp method:
EDR_ramp=Rpeak_EDR(ECG,R,fs);
% K-PCA method:
% EDR_KPCA=KPCA_EDR(ECG,R,fs);
% Respiratory sinus arrhythmia (RSA)
% EDR_RSA=RSA_resp(ECG,R,fs);


% RESULTS:
plot(EDR_physionet(:,1),EDR_physionet(:,2)),hold on, 
plot(EDR_ramp(:,1),EDR_ramp(:,2),'k'),
% plot(EDR_KPCA(:,1),100*EDR_KPCA(:,2),'r'),
% plot(EDR_RSA(:,1),100*EDR_RSA(:,2),'m'),xlabel('sec')
legend('physionet EDR','EDR-RAMP','EDR-KPCA','RSA'),