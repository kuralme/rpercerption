function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

n = length(x1);
A = zeros(n,9);

for m=1:n
    A(m,:)=[[x1(m,1)*x2(m,1)],...
            [x1(m,1)*x2(m,2)],...
            [x1(m,1)],...
            [x1(m,2)*x2(m,1)],...
            [x1(m,2)*x2(m,2)],...
            [x1(m,2)],...
            [x2(m,1)^2],...
            [x2(m,2)^2],...
            [1]];    
end

[~,~,V] = svd(A);
H = reshape(V(:,end),3,3)';

% SVD cleanup
if (rank(H)==3)
  [U,D,V] = svd(H);
  D(3,3) = 0;
  F2 = U * D * V.';
  F = F2 / norm(F2);
else
  F = H;
end
end

  



