function [O] = local_frame(u, v, rot_y2x)
% construct a normal orthogonal frame from a local element
% in global frame
% Inputs:
% u, v: 1x3 vector in global xyz frame
% rot_y2x: optional rotation angle
assert(length(u) == 3 && length(v) == 3);
if isempty(rot_y2x)
    rot_y2x = 0;
end

L = norm(v-u);
c = (v-u) / L;
O = zeros(3,3);

if abs(c(3))==1.0
    O(3,1) = 1.0;
    O(1:2,2:3) = [cos(rot_y2x), -sin(rot_y2x); sin(rot_y2x), cos(rot_y2x)];
else
    new_x = c;
    new_y = cross([0,0,1], new_x);
    R = rot_axis(new_x, rot_y2x);
    
    O(1,:) = new_x;
    O(2,:) = new_y;
    O(3,:) = cross(new_x, new_y);
end

end