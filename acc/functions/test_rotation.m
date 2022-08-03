clear all; close all; clc;

p0 = [0.0566, 0.0232, 1.0264]';
p1 = [0.0156, 1.0496, -0.0037]'; 

p0 = p0/norm(p0);
% p1_noise = p1 + normrnd(0, 0.3, [3, 1]);
p1 = p1/norm(p1);
% p1_noise = p1_noise/norm(p1_noise);

GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;
              norm(cross(A,B)) dot(A,B)  0;
              0              0           1];

FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

UU = @(Fi,G) Fi*G*inv(Fi);

rot = UU(FFi(p0, p1), GG(p0, p1));
% rot_noise = UU(FFi(p0, p1_noise), GG(p0, p1_noise));


yaw = atan2d(rot(2, 1), rot(1, 1));
pitch = atan2d(-rot(3, 1), sqrt(rot(3, 2)^2 + rot(3, 3)^2));
roll = -atan2d(rot(3, 2), rot(3, 3));