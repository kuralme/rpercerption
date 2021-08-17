function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

% ^ Means, use xc = inv(K)*x ^

% Making the world points and image points homogenous
X(:, end+1) = 1;  % Nx4
x(:, end+1) = 1;  % Nx3

for i=1:length(X)
    A(i*2-1,1:4)= zeros(1,4);
    A(i*2-1,5:8)= -X(i,:);
    A(i*2-1,9:12)= x(i,2)*X(i,:);
    
    A(i*2,1:4) = X(i,:);
    A(i*2,5:8) = zeros(1,4);
    A(i*2,9:12)= -x(i,1)*X(i,:);
end

% Ax = 0
[~,~,Vn] = svd(A);
P = reshape(Vn(:,12)/Vn(12,12),4,3).';

%R = (K^-1)*P
P = inv(K)*P; %calibrated
[U,D,V] = svd(P(:,1:3));

% R and t with cleanup
if ( det(U*V.')>0 )
    R = U*V.'; % R+
    t = P(:,4)/D(1,1);
    C = -R.'*t;
elseif( det(U*V.')<0 )
    R = -U*V.'; % R+
    t = -P(:,4)/D(1,1);
    C = -R.'*t;
end
end    





