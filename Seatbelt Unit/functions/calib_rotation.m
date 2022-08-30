function [R, t] = calib_rotation(main, ref)
% main: 3xN matrix
% ref: 3xN matrix
[~, n] = size(main);
centroid_main = mean(main, 2);
centroid_ref = mean(ref, 2);

main_m = main - repmat(centroid_main, 1, n);
ref_m = ref - repmat(centroid_ref, 1, n);

H = main_m*ref_m';

[U, ~, V] = svd(H);
R = V*U';
if det(R) < 0
    V(:, 3) = V(:, 3)*-1;
    R = V*U';
end

t = -R*centroid_main + centroid_ref;



% [~, n] = size(main);
% 
% rot_avg = zeros(3, 3);
% for i = 1:n
%     main_single = main(:, i);
%     ref_single = ref(:, i);
%     main_single = main_single/norm(main_single);
%     ref_single = ref_single/norm(ref_single);
%     
%     GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;
%                   norm(cross(A,B)) dot(A,B)  0;
%                   0              0           1];
%     
%     FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
%     
%     UU = @(Fi,G) Fi*G*inv(Fi);
%     
%     rot = UU(FFi(main_single, ref_single), GG(main_single, ref_single));
%     rot_avg = rot_avg + rot;
% end
% 
% rot_avg = rot_avg/n;
end