function Rn=KPCA_EDR(ECG,ECGR,fs)
% ECG-derived Respiratory (EDR) signal from a single-lead ECG signal, using
% Kernel-PCA method;
%
% Inputs:
%
% ECG
%       A 1xN integer array containing the ECG signal 
%
% ECGR
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
%
% written based on:
%
% Widjaja, Devy, et al. "Application of kernel principal component analysis 
% for single-lead-ECG-derived respiration." Biomedical Engineering, 
% IEEE Transactions on 59.4 (2012): 1169-1176.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. 
% 
%--------------------------------------------------------------------------
%% Initialization 
W=[.1,.45];% respiration frequency band
m=121; % window length
% fixing the orientation of vectors
ECGR=round(ECGR*fs);
if size(ECG,1)<size(ECG,2),
    ECG=ECG';
end
if size(ECGR,1)<size(ECGR,2),
    ECGR=ECGR';
end
n=length(ECGR);
%%%%%%%%%%%%%%%%
%% Construction of Input matrix (X)
X=zeros(m,n);
for i=1:n
    if ECGR(i)>=(1+(m-1)/2) && ECGR(i)+((m-1)/2)<=length(ECG)
    X(:,i)=ECG(ECGR(i)-((m-1)/2):ECGR(i)+((m-1)/2),1);
    elseif ECGR(i)<(1+(m-1)/2)
        X(1:(1+(m-1)/2)-ECGR(i),i)=ECG(i)*ones((1+(m-1)/2)-ECGR(i),1);
        X((2+(m-1)/2)-ECGR(i):m,i)=ECG(1:ECGR(i)+((m-1)/2),1);
    elseif ECGR(i)+((m-1)/2)>length(ECG)
        X((length(ECG)-ECGR(i)+(2+(m-1)/2)):m,i)=ECG(length(ECG))*ones(m-(length(ECG)-ECGR(i)+(1+(m-1)/2)),1);
        X(1:(length(ECG)-ECGR(i)+(1+(m-1)/2)),i)=ECG(ECGR(i)-((m-1)/2):length(ECG),1);
    end
end
%%%%%%%%%%%%%%%%
%% KPCA and choosing the best sigma for Kernel
reso=1e7; % steps for sigma^2  
s1=m*mean(var(X));s=s1/100:reso:s1*100; % a range of sigma^2 to be tested

for i=1:length(s)
   [eigval(:,i), eigvec{i}] = kpca(X, 'RBF_kernel', s(i) ,[],'svd',10,'o'); 
   edif(i)=eigval(1,i)-sum(eigval(2:m,i)); % to be used later to optimize sigma^2
end
[mi,ii]=max(edif); % ii is the index for optimal sigma^2, i.e. signa^2(ii)
EigVal=eigval(:,ii);EigVec=eigvec{ii};
%%%%%%%%%%%%%%%
%% Reconstruction
size(X)
Xr=preimage_rbf(X,s(ii),EigVec,X,'r',10,1000);
[mxr1,mxr2]=min(abs(mean(Xr')));
R=Xr(mxr2,:);
Rn(:,2)=-(R-mean(R))';Rn(:,1)=ECGR'/fs;
%R=Xr(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% respiration rate
% B1=interp1(ECGR,R,1:1000/4:length(ECG),'cubic');
% lB=length(B1);
% 
% f = [0 .9*W(1) W(1) W(2) 1.1*W(2) 2]/2;
% a = [0 0 1 1 0 0] ;
% b = firls(30,f,a);
% B1=filtfilt(b,1,B1);
% 
% bb=findpeaks1(B1);
% b1=B1(bb);
% B1=B1/mean(b1(b1>prctile(b1,75)));
% B(1,:)=(1:1000/4:length(ECG))/1000;
% B(2,:)=B1;
% ra(1,:)=bb(B1(bb)>0.3);
% Rra=diff(ra);
% rxa=ra;
% clear ra
% ra=rxa(1:length(Rra));
% ra(2,:)=Rra;
% % % to 1HZ
% % ra(2,:)=4*ra(2,:);
% % for i=1:floor(((lB-80)/20)+1)
% % a(i,:) = arcov(B1(20*(i-1)+1:20*(i-1)+80),40);
% % a(i,[1:3,38:41])=0;
% % c(i)=.95*max(a(i,:));
% % pa=findpeaks1(a(i,:));
% % ra(i)=min(pa(a(i,pa)>=c(i)));
% % ra(i)=(ra(i)-1)/40;
% % clear pa
% % end
%  Ra=mean(ra(2,:));