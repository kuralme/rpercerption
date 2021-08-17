function E = EssentialMatrixFromFundamentalMatrix(F,K)
%% EssentialMatrixFromFundamentalMatrix
% Use the camera calibration matrix to esimate the Essential matrix
% Inputs:
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     F - size (3 x 3) fundamental matrix from EstimateFundamentalMatrix
% Outputs:
%     E - size (3 x 3) Essential matrix with singular values (1,1,0)


E = K.'*F*K;

[U,D,V] = svd(E);
if (D(1,1)~= D(2,2))
  D = diag([1,1,0]);
end

E2 = U * D * V.';
E = E2 / norm(E2);

end

