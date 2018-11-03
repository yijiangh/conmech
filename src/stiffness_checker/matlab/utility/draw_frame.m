function [] = draw_frame(N, T, S, F, D, draw_force, thkMin, thkMax, magn)
% Input
%
% N = node coordinates
% (number of nodes)-by-dim matrix with N(n,1) and N(n,2) the X and
% Y coordinates of node n
%
% T = truss topology
% (number of elements)-by-2 matrix with T(e,1) and T(e,2) the indices of
% the starting and ending nodes of element e
%
% F = element forces
% (number of elements)-by-1 matrix with F(e,1) the force in element e.
% May be empty.
%
% thkMin = line thickness used to draw elements with zero force.
% May be empty.
%
% thkMax = line thickness used to draw elements with maximum force.
% May be empty.
%
% D = node displacements
% (number of nodes)-by-2 matrix with D(n,1) and D(n,2) the displacements of
% node n in the X and Y directions.
% May be empty.
%
% magn = displacement magnification factor.
% May be empty.

dim = size(N,2);
nNodes = size(N,1);
nElements = size(T,1);

% Verify and complete line thickness input
if isempty(thkMin) == 0
    sizeOfThk = size(thkMin);
    if sizeOfThk(1,1) ~= 1 || sizeOfThk(1,2) ~= 1
        error('Line thickness input must be empty or a number')
    end
end
if isempty(thkMax) == 0
    sizeOfThk = size(thkMax);
    if sizeOfThk(1,1) ~= 1 || sizeOfThk(1,2) ~= 1
        error('Line thickness input must be empty or a number')
    end
end
if isempty(thkMin) && isempty(thkMax)
    thkMin = 1;
    thkMax = 1;
elseif isempty (thkMax)
    thkMax = thkMin;
elseif isempty (thkMin)
    thkMin = thkMax;
elseif thkMin > thkMax
    temp = thkMin;
    thkMin = thkMax;
    thkMax = temp;
end

% Determine line thicknesses and colors
if isempty(draw_force) || 0 == draw_force || isempty(F)
    thicknesses = zeros(1,nElements);
    colors = zeros(nElements,3);
    if isempty(F)
        thk = 0.5*(thkMin+thkMax);
        for e=1:1:nElements
            thicknesses(1,e) = thk;
            colors(e,:) = [0,0,0]; % black
        end
    else
        sizeOfF = size(F);
        if sizeOfF(1,1) ~= nElements || sizeOfF(1,2) ~= 1
            error('Invalid element force input')
        end
        
        thkFact = (thkMax-thkMin)/max(max(F), -min(F));
        for e=1:1:nElements
            thicknesses(1,e) = thkMin + thkFact*abs(F(e,1));
            if (F(e,1) < 0)
                colors(e,:) = [1,0,0]; % red
            else
                colors(e,:) = [0,0,1]; % green
            end
        end
    end
end

% Verify and complete node displacements input
if isempty(D) == 0
    sizeOfD = size(D);
    if sizeOfD(1,1) ~= nNodes || sizeOfD(1,2) ~= 2
        error('Invalid node displacement input')
    end
    
    if isempty(magn)
        magn = 1;
    else
        sizeOfMagn = size(magn);
        if sizeOfMagn(1,1) ~= 1 || sizeOfMagn(1,2) ~= 1
            error('Displacement magnification input must be empty or a number')
        end
    end
    
    for n=1:1:nNodes
        N(n,1) = N(n,1)+magn*D(n,1);
        N(n,2) = N(n,2)+magn*D(n,2);
    end
end

% Set axes
h = figure(1);
minX = min(N(:,1));
maxX = max(N(:,1));
minY = min(N(:,2));
maxY = max(N(:,2));

if 3 == dim
    minZ = min(N(:,3));
    maxZ = max(N(:,3));
    gap = 0.1*max([maxX-minX, maxY-minY, maxZ-minZ]);
    
    axis equal
    axis([minX-gap, maxX+gap, minY-gap, maxY+gap, minZ-gap, maxZ+gap])
    xlabel('x axis');
    ylabel('y axis');
    zlabel('z axis');
    view(26,10.35);
else
    % dim = 2
    axis equal
    axis([minX-gap, maxX+gap, minY-gap, maxY+gap])
    xlabel('x axis');
    ylabel('y axis');
end

hold on;
alpha = 0.2;
plot_undeformed_frame(N,T,dim);
plot_fixities(N,S,dim,alpha);

hold off;

end

function plot_undeformed_frame(N, T, dim)
nElements = size(T,1);
undeformed_color = [0.25, 0.25, 0.25];
for e=1:1:nElements
    if 3 == dim
        plot3(...
            [N(T(e,1),1);N(T(e,2),1)],...
            [N(T(e,1),2);N(T(e,2),2)],...
            [N(T(e,1),3);N(T(e,2),3)],...
            'Color',undeformed_color,'LineWidth',1.2);
    else
        line(...
            [N(T(e,1),1);N(T(e,2),1)],...
            [N(T(e,1),2);N(T(e,2),2)],...
            'Color',undeformed_color,'LineWidth',1.2);
    end
end
end

function plot_deformed_frame()

end

function plot_fixities(N, S, dim, alpha)
n_Fix = size(S,1);
xp = N(S(:,1),:);
fix_color = [0.4660, 0.6740, 0.1880];
lw = 1.2;
c_alpha = 0.4;
cir_a = 0.5;

if 3 == dim
    for s=1:1:n_Fix
        p = xp(s,:);
        scatter3(p(1),p(2),p(3),'filled','k');
        TR = S(s,2:7);
        
        if TR(1)
            plot3(...
                [p(1)-(alpha); p(1)],...
                [p(2); p(2)],...
                [p(3); p(3)],'Color',fix_color,'LineWidth',lw);
        end
        
        if TR(2)
            plot3(...
                [p(1); p(1)],...
                [p(2)-(alpha); p(2)],...
                [p(3); p(3)],'Color',fix_color,'LineWidth',lw);
        end
        
        if TR(3)
            plot3(...
                [p(1); p(1)],...
                [p(2); p(2)],...
                [p(3)-(alpha); p(3)],'Color',fix_color,'LineWidth',lw);
        end
        
        if TR(4)
            plot_circle(p-cir_a*alpha*[1,0,0], [1,0,0], alpha*c_alpha, 3, fix_color, lw);
        end
        if TR(5)
            plot_circle(p-cir_a*alpha*[0,1,0], [0,1,0], alpha*c_alpha, 3, fix_color, lw);
        end
        if TR(6)
            plot_circle(p-cir_a*alpha*[0,0,1], [0,0,1], alpha*c_alpha, 3, fix_color, lw);
        end
    end
else
    % dim = 2
    for s=1:1:n_Fix
        p = xp(s,:);
        scatter(p(1),p(2),'filled','k');
        TR = S(s,2:7);
        
        if TR(1)
            line(...
                [p(1)-(alpha); p(1)],...
                [p(2); p(2)],'Color',fix_color,'LineWidth',lw);
        end
        
        if TR(2)
            line(...
                [p(1); p(1)],...
                [p(2)-(alpha); p(2)],'Color',fix_color,'LineWidth',lw);
        end
        
        if TR(6)
            plot_circle(p, [0,0,1], alpha*c_alpha, 2, fix_color, lw);
        end
    end
end

end

function plot_load()

end