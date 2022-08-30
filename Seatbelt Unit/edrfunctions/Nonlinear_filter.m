function data = Nonlinear_filter(x)

% //Copyright (c) 2012, Chengyu Liu
% //All rights reserved.

% $ This function is used to reduce the gauss noise existing in a physioilogical signal using the nonlinear filter.
%   Input x: raw signal
%   Output data: filtered signal

% $ Author:  Chengyu Liu (bestlcy@sdu.edu.cn) 
%            Institute of Biomedical Engineering,
%            Shandong University
% $ Last updated:  2012.08.02

x          = x-mean(x);
tt         = 0.0001:0.0001:1;
indx       = x;
indx(x>0)  = 1;
indx(x<=0) = -1;
mm         = max(abs(x));
x_tt       = abs(x)/mm;

%% generate the nonlinear filter using 3 segment functions
M            = 2.2;
ratio        = 3;
t1           = 0.0001:0.0001:1;
y1           = t1.^3;
t2           = 0.0001:0.0001:(M-2);
y2           = t2*ratio+1;
y3           = y1;
y3(end:-1:1) = y3;
y3           = -y3;
y3           = y3+max(y2)+1;
y            = [y1 y2 y3];
yy           = resample(y,11000,22000);
yy           = yy(1:10000);
Trans_fun2   = yy/max(yy);

%% nonlinear transfering
for i = 1:length(x_tt)
    [indd indx_cur] = min(abs(x_tt(i)-tt));
    x(i)            = Trans_fun2(indx_cur)*indx(i);
end
data = x*mm;
