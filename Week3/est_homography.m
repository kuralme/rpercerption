function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE% 

% A= [
% -x1  -y1  -1   0    0    0   x1*xp1   y1*xp1   xp1;
%  0    0    0 -x1   -y1  -1   x1*yp1   y1*yp1   yp1;
% -x2  -y2  -1   0    0    0   x2*xp2   y2*xp2   xp2;
%  0    0    0 -x2   -y2  -1   x2*yp2   y2*yp2   yp2;
% -x3  -y3  -1   0    0    0   x3*xp3   y3*xp3   xp3;
%  0    0    0 -x3   -y3  -1   x3*yp3   y3*yp3   yp3;
% -x4  -y4   -1  0    0    0   x4*xp4   y4*xp4   xp4;
%  0    0    0  -x4  -y4  -1   x4*yp4   y4*yp4   yp4];
%

%ax: 4x9
ax =[[-video_pts(:,1)],...
     [-video_pts(:,2)],...
     [-1*ones(4,1)],...
     [0*ones(4,1)],...
     [0*ones(4,1)],...
     [0*ones(4,1)],...
     [video_pts(:,1).*logo_pts(:,1)],...
     [video_pts(:,2).*logo_pts(:,1)],...
     [logo_pts(:,1)]];

%ay: 4x9   
ay =[[0*ones(4,1)],...
     [0*ones(4,1)],...
     [0*ones(4,1)],...
     [-video_pts(:,1)],...
     [-video_pts(:,2)],...
     [-1*ones(4,1)],...
     [video_pts(:,1).*logo_pts(:,2)],...
     [video_pts(:,2).*logo_pts(:,2)],...
     [logo_pts(:,2)]];

%A: 8x9
A = [ax(1,:);
     ay(1,:);
     ax(2,:);
     ay(2,:);
     ax(3,:);
     ay(3,:);
     ax(4,:);
     ay(4,:)];
 
% ****Another method for A*****
% h = zeros(9,1);
% A = zeros(8,9);
% 
% for m=1:4
%    A(2*m-1,1)= -video_pts(m,1);
%    A(2*m-1,2)= -video_pts(m,2);
%    A(2*m-1,3)= -1;
%    A(2*m-1,4)= 0;
%    A(2*m-1,5)= 0;
%    A(2*m-1,6)= 0;
%    A(2*m-1,7)= video_pts(m,1) * logo_pts(m,1);
%    A(2*m-1,8)= video_pts(m,2) * logo_pts(m,1);
%    A(2*m-1,9)= logo_pts(m,1);
%    A(2*m,1)= 0;
%    A(2*m,2)= 0;
%    A(2*m,3)= 0;
%    A(2*m,4)= -video_pts(m,1);
%    A(2*m,5)= -video_pts(m,2);
%    A(2*m,6)= -1;
%    A(2*m,7)= video_pts(m,1) * logo_pts(m,2);
%    A(2*m,8)= video_pts(m,2) * logo_pts(m,2);
%    A(2*m,9)= logo_pts(m,2);
% end
% *************************

[~,~,V] = svd(A);

H = reshape(V(:,end),3,3).';

end

