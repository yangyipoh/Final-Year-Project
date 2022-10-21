function [R, t] = calib_rotation(A, B)
% main: 3xN matrix
% ref: 3xN matrix
[num_rows, num_cols] = size(A);
if num_rows ~= 3
    error("matrix A is not 3xN, it is %dx%d", num_rows, num_cols)
end

[num_rows, num_cols] = size(B);
if num_rows ~= 3
    error("matrix B is not 3xN, it is %dx%d", num_rows, num_cols)
end

    % find mean column wise
centroid_A = mean(A, 2);
centroid_B = mean(B, 2);

    % subtract mean
Am = A - repmat(centroid_A, 1, num_cols);
Bm = B - repmat(centroid_B, 1, num_cols);

    % calculate covariance matrix (is this the correct terminology?)
H = Am * Bm';

    %if rank(H) < 3
    %    error(sprintf("rank of H = %d, expecting 3", rank(H)))
    %end

% find rotation
[U,S,V] = svd(H);
R = V*U';

if det(R) < 0
%     fprintf("det(R) < R, reflection detected!, correcting for it ...\n");
    V(:,3) = V(:,3) * -1;
    R = V*U';
end

t = -R*centroid_A + centroid_B;
end