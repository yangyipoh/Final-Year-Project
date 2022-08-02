clear all; close all; clc;

p0 = [4, -2, 1]';
p0 = p0/norm(p0);
p1 = [1, -1, 3]'; 
p1_noise = p1 + normrnd(0, 0.3, [3, 1]);
p1 = p1/norm(p1);
p1_noise = p1_noise/norm(p1_noise);

GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;
              norm(cross(A,B)) dot(A,B)  0;
              0              0           1];

FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

UU = @(Fi,G) Fi*G*inv(Fi);

rot = UU(FFi(p0, p1), GG(p0, p1));
rot_noise = UU(FFi(p0, p1_noise), GG(p0, p1_noise));
