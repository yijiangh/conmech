function [R] = local_frame(u, v, rot_y2x)
% construct a normal orthogonal frame from a local element
% in global frame
% Inputs:
% u, v: 1x(2 or 3) vector in global xyz frame
% rot_y2x: optional rotation angle

assert(length(u) == length(v) && size(u,2) == size(v,2));
dim = size(u,2);

assert(~all(u == v));

if isempty(rot_y2x)
    rot_y2x = 0;
end

L = norm(v-u);
c = (v-u) / L;
Cx = c(1);
Cy = c(2);

if 3 == dim
    Cz = c(3);
    R = zeros(3,3);
    
    if abs(Cz)==1.0        
        R(1,3) = -Cz;
        R(2,2) = 1;
        R(3,1) = Cz;
        
        R = R * rot_axis([1,0,0], rot_y2x);
        R = R';
    else
        new_x = c;
        
        new_y = cross(new_x, [0,0,1]);
        new_y = -new_y/norm(new_y);
        
        R(:,1) = new_x';
        R(:,2) = new_y';
        R(:,3) = cross(new_x, new_y)';
        
        R = R * rot_axis([1,0,0], rot_y2x);
        
        % THIS IS ESSENTIAL!!!
        R = R';
    end
else
    R = zeros(3,3);
    R(1,1) = Cx;
    R(1,2) = Cy;
    R(2,1) = -Cy;
    R(2,2) = Cx;
    R(3,3) = 1;
end

end