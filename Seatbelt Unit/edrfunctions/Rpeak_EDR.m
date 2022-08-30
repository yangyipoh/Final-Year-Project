function Rn=Rpeak_EDR(ECG,ECGR,fs)

% ECG-derived Respiratory (EDR) signal from a single-lead ECG signal, based 
% on the R-peaks and using RAMP method;
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

if size(ECG,1)>size(ECG,2);
    ECG=ECG';
end
if size(ECGR,1)>size(ECGR,2);
    ECGR=ECGR';
end

ECGR=round(ECGR*fs);

y = medfilt1(ECG,200);
Y=medfilt1(y,600);
ECGb=ECG-Y;
for i=1:length(ECGR)
R(i)=trapz(ECGb(ECGR(i)+1:min(ECGR(i)+100,length(ECGb))));
end
Rn(:,2)=-(R-mean(R))';Rn(:,1)=ECGR'/fs;
