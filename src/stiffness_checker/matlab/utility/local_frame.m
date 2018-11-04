function [R] = local_frame(u, v, rot_y2x)
% construct a normal orthogonal frame from a local element
% in global frame
% Inputs:
% u, v: 1x(2 or 3) vector in global xyz frame
% rot_y2x: optional rotation angle

assert(length(u) == length(v) && size(u,2) == size(v,2));
dim = size(u,2);

if isempty(rot_y2x)
    rot_y2x = 0;
end

L = norm(v-u);
c = (v-u) / L;

if 3 == dim
    R = zeros(3,3);
    
    if abs(c(3))==1.0
        R(3,1) = 1.0;
        R(1:2,2:3) = [cos(rot_y2x), -sin(rot_y2x); sin(rot_y2x), cos(rot_y2x)];
    else
        new_x = c;
        new_y = cross([0,0,1], new_x);
        R = rot_axis(new_x, rot_y2x);
        
        R(1,:) = new_x;
        R(2,:) = new_y;
        R(3,:) = cross(new_x, new_y);
    end
else
    R = zeros(3,3);
    R(1,1) = c(1);
    R(1,2) = c(2);
    R(2,1) = -c(2);
    R(2,2) = c(1);
    R(3,3) = 1;
end

end