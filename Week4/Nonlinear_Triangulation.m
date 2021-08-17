function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 

for i =1:length(X0)
    re = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), X0(i,:).');  
    re2 = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), re); 
    re3 = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), re2); 
    X(i,:) = re3.';
end
end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)

% Estimated homogenous image points
x1est = K*R1*(X0-C1);
x2est = K*R2*(X0-C2);
x3est = K*R3*(X0-C3);

b = [ x1(1,1);
      x1(1,2);
      x2(1,1);
      x2(1,2);
      x3(1,1);
      x3(1,2)];

fX = [ x1est(1)/x1est(3);
       x1est(2)/x1est(3);
       x2est(1)/x2est(3);
       x2est(2)/x2est(3);
       x3est(1)/x3est(3);
       x3est(2)/x3est(3)];
   
J = [Jacobian_Triangulation(C1, R1, K, X0);
     Jacobian_Triangulation(C2, R2, K, X0);
     Jacobian_Triangulation(C3, R3, K, X0)];
 
% Update: X = X + deltaX
X = X0 + inv(J.'*J)*J.'*(b-fX);
end

function grad_fX = Jacobian_Triangulation(C, R, K, X)
x  = K*R*(X-C);
f  = K(1,1);
px = K(1,3);
py = K(2,3);

grad_uX = [f*R(1,1)+px*R(3,1),  f*R(1,2)+px*R(3,2),  f*R(1,3)+px*R(3,3)];
grad_vX = [f*R(2,1)+py*R(3,1),  f*R(2,2)+py*R(3,2),  f*R(2,3)+px*R(3,3)];
grad_wX = R(3,:);

grad_fX = [ (x(3)*grad_uX-x(1)*grad_wX)/(x(3)*x(3));
            (x(3)*grad_vX-x(2)*grad_wX)/(x(3)*x(3))];
end
