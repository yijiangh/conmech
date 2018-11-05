function plot_circle(center,z_axis,r,dim,color,lw)
% Input:
%   center: 1x3 vector for circle center
%   z_axis: 1x3 vector
%   r: circle radius
%   dim: dimension (2 or 3)
%   color:
%   lw: line width

if 2 == dim
    theta=0:0.1:2*pi+0.1;
    x=r*cos(theta)+center(1);
    y=r*sin(theta)+center(2);
    
    plot(x,y,'Color',color,'LineWidth',lw);
else
    O = local_frame(center, center+z_axis, 0);
    p = r*O(:,2);
    
    xp=[];
    for theta=0:0.1:2*pi+0.1
        xp = [xp; center + (rot_axis(z_axis, theta)*p)'];
    end
    plot3(...
        xp(:,1),...
        xp(:,2),...
        xp(:,3),'Color',color,'LineWidth',lw);
end