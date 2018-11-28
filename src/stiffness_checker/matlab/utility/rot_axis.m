function [R] = rot_axis(r_v, theta)
% return 3x3 rotation matrix for rotation around
% 3d axis
% Input:
%   r_v: rotation axis vector
%   theta: rotation angle in radian
% see: http://ksuweb.kennesaw.edu/~plaval/math4490/rotgen.pdf

r_v = r_v / norm(r_v);
R = zeros(3,3);
C = cos(theta);
S = sin(theta);
if theta == pi
    C=round(C);
    S=round(S);
end
t = 1-C;

ux = r_v(1);
uy = r_v(2);
uz = r_v(3);

R(1,1) = t*ux^2 + C;
R(1,2) = t*ux*uy - S*uz;
R(1,3) = t*ux*uz + S*uy;

R(2,1) = t*ux*uy + S*uz;
R(2,2) = t*uy^2 + C;
R(2,3) = t*uy*uz - S*ux;

R(3,1) = t*ux*uz - S*uy;
R(3,2) = t*uy*uz + S*ux;
R(3,3) = t*uz^2 + C;

end