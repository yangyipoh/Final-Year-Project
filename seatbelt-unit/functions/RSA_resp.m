function Rn=RSA_resp(ECG,mR,fs)
% ECG-derived Respiratory basde on Respiratory sinus arrhythmia (RSA);
%
% Inputs:
%
% ECG
%       A 1xN integer array containing the ECG signal 
%
% mR
%       A 1xM vctor containing the locations of r peaks on signal in
%       sec.
% fs
%       sampling frequency
%
% output:
%
% Rn
%       A Mx2 matrix containing time in sec and edr
%
% -----------------------------------------------------------------------
% Released under the GNU General Public License
%
% Copyright (C) 2014 Faezeh Marzbanrad
% The University of Melbourne
% 
% faezehm@student.unimelb.edu.au
% fmrad64@gmail.com 

if size(ECG,1)>size(ECG,2);
    ECG=ECG';
end
if size(mR,1)>size(mR,2);
    mR=mR';
end


W=[.1,.45];
mR=round(mR*fs);
mRR=diff(mR);
%mRR=Adaptivepreproc(mRR,0.5);
% Interpolation 10 HZ
B1=interp1(mR(1:length(mRR)),mRR,1:1000/10:length(ECG),'cubic');
B1 = fillmissing(B1, 'linear');
% band-pass filtered 0.1–0.45 Hz by a least-square FIR filter (fall-off:60 dB)
f = [0 .9*W(1) W(1) W(2) 1.1*W(2) 5]/5;
a = [0 0 1 1 0 0] ;
b = firls(30,f,a);
B=filtfilt(b,1,B1);
% Standardization by dividing the time series by the 75%-percentil of all local maxima.
bb=findpeaks1(B);
b1=B(bb);
R=(1:1000/10:length(ECG))/1000;
R(2,:)=B/mean(b1(b1>prctile(b1,75)));
Rn=R';
Rn(:,2)=Rn(:,2)-mean(Rn(:,2));
% ra(1,:)=bb(B(bb)>0.3);
% ra=ra(1:(length(ra)-1));
% ra(2,:)=diff(bb(B(bb)>0.3));
% % to 1HZ scale
% ra(2,:)=10*ra(2,:);
% Ra=mean(ra(2,:));