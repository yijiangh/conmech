function [] = draw_frame(N, T, S, Load, F, R, D, thkMin, thkMax, magn, r_scale)
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
% S = fixities
% (number of fixities)-by-(4 or 7) matrix.
%
% Load = load
% (number of fixities)-by-(4 or 7) matrix.
%
% F = element forces
% (number of elements)-by-(4 or 7) matrix with F(e,1) the force in element e.
% May be empty.
% In 2D, this includes element's axial force, shear force, and moment (3)
% In 3D, this includes element's axial force, shear force (Vy, Vz),
% moment (My, Mz), and torsion (Tx).
%
% D = node displacements
% (number of nodes)-by-dim matrix with D(n,1) and D(n,2) the displacements of
% node n in the X and Y directions.
% May be empty.
%
% thkMin = line thickness used to draw elements with zero force.
% May be empty.
%
% thkMax = line thickness used to draw elements with maximum force.
% May be empty.
%
% magn = displacement magnification factor.
% May be empty.

dim = size(N,2);
if 2 == dim
    assert(size(S,2) == (1+3));
    %     assert(isempty(F) || size(F,2) == (1+3));
else
    assert(size(S,2) == (1+6));
    %     assert(isempty(F) || size(F,2) == (1+6));
end

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
undeformed_color = [0.25, 0.25, 0.25];
thicknesses = zeros(1,nElements);
colors = zeros(nElements,3);
if isempty(F)
    thk = 0.5*(thkMin+thkMax);
    for e=1:1:nElements
        thicknesses(1,e) = thk;
        colors(e,:) = undeformed_color; % black
    end
else
    sizeOfF = size(F);
    if sizeOfF(1,1) ~= nElements
        error('Invalid element force input')
    end
    
    thkFact = (thkMax-thkMin)/max(max(F(:,1)), -min(F(:,1)));
    for e=1:1:nElements
        thicknesses(1,e) = thkMin + thkFact*abs(F(e,1));
        if dim == 2
            if (F(e,1) > 0)
                colors(e,:) = [1,0,0]; % red
            else
                colors(e,:) = [0,0,1]; % blue
            end
        else
            if (F(e,1) > 0)
                colors(e,:) = [1,0,0]; % red
            else
                colors(e,:) = [0,0,1]; % blue
            end
        end
    end
end

% Set axes
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
    view(10,10.35);
else
    % dim = 2
    gap = 0.1*max([maxX-minX, maxY-minY]);
    axis equal
    axis([minX-gap, maxX+gap, minY-gap, maxY+gap])
    xlabel('x axis');
    ylabel('y axis');
end

if dim == 2
    N = N(:,1:2);
else
    N = N(:,1:3);
    D = D(:,1:3);
end

h = figure(1);
hold on;
alpha = 1;
plot_frame(N,[],T,undeformed_color,[],dim);
% plot_fixities(N,S,dim,alpha);

if ~isempty(Load)
    plot_load(N,Load,dim,alpha*r_scale)
end

if ~isempty(R)
    plot_reaction(N,R,dim,alpha*r_scale)
end

plot_frame(N,D,T,colors,thicknesses,dim,magn);

hold off;

end

function plot_frame(N, D, T, colors, thickness, dim,magn)
nElements = size(T,1);
if 1 == size(colors,1)
    c = colors;
    for e=1:1:size(T,1)
        colors(e,:) = c;
    end
end
if isempty(thickness) || 1 == length(thickness)
    thickness = ones(nElements)*1.2;
end
if ~isempty(D)
    if size(D,1)~=size(N,1)
        error('Invalid node displacement input')
    end
end

for e=1:1:nElements
    if 3 == dim
        plot3(...
            [N(T(e,1),1);N(T(e,2),1)],...
            [N(T(e,1),2);N(T(e,2),2)],...
            [N(T(e,1),3);N(T(e,2),3)],...
            'Color',colors(e,:),'LineWidth',thickness(e));
    else
        if ~isempty(D)
            if size(D,2)>2
                % contains rotation, beam or frame
                D_beam = draw_cubic_bent_beam(N(T(e,1),:), N(T(e,2),:), ...
                    D(T(e,1),:), D(T(e,2),:), magn);

                for b_i=1:1:size(D_beam,1)-1
                    line([D_beam(b_i,1); D_beam(b_i+1,1)],...
                        [D_beam(b_i,2); D_beam(b_i+1,2)],...
                        'Color',colors(e,:),'LineWidth',thickness(e));
                end
            else
                % truss, only translation, no beam shape intepolation
                if isempty(magn)
                    magn = 1;
                end
                
                for n=1:1:size(N,1)
                    N(n,:) = N(n,:) + magn*D(n,:);
                end
                line(...
                    [N(T(e,1),1);N(T(e,2),1)],...
                    [N(T(e,1),2);N(T(e,2),2)],...
                    'Color',colors(e,:),'LineWidth',thickness(e));
            end
        else
            line(...
                [N(T(e,1),1);N(T(e,2),1)],...
                [N(T(e,1),2);N(T(e,2),2)],...
                'Color',colors(e,:),'LineWidth',thickness(e));
        end
    end
end
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
            plot_circle(p-cir_a*alpha*[1,0,0], [1,0,0], alpha*c_alpha, dim, fix_color, lw);
        end
        if TR(5)
            plot_circle(p-cir_a*alpha*[0,1,0], [0,1,0], alpha*c_alpha, dim, fix_color, lw);
        end
        if TR(6)
            plot_circle(p-cir_a*alpha*[0,0,1], [0,0,1], alpha*c_alpha, dim, fix_color, lw);
        end
    end
else
    % dim = 2
    for s=1:1:n_Fix
        p = xp(s,:);
        scatter(p(1),p(2),'filled','k');
        TR = S(s,2:4);
        
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
        
        if TR(3)
            plot_circle(p, [], alpha*c_alpha, dim, fix_color, lw);
        end
    end
end

end

function plot_load(N, Load, dim, alpha)
assert(any(Load(:,1)) <= size(N,1));
load_color = [0.3010, 0.7450, 0.9330];
lw = 2;
n_Load = size(Load, 1);

for f=1:1:n_Load
    switch dim
        case 2
            quiver(N(Load(f,1),1),N(Load(f,1),2),...
                alpha*Load(f,2),alpha*Load(f,3),...
                'Color', load_color, 'LineWidth', lw);
            if size(Load,2)-1 == 3 && Load(f,4) ~= 0
                plot_circle(N(Load(f,1),:), [], alpha*0.5, dim, load_color, lw);
            end
        case 3
            quiver3(N(Load(f,1),1),N(Load(f,1),2),N(Load(f,1),3),...
                alpha*Load(f,2),alpha*Load(f,3),alpha*Load(f,4),...
                'Color', load_color, 'LineWidth', lw);
            if size(Load,2)-1 == 6
                E = eye(3);
                for s=1:1:3
                    if Load(f,4+s) ~= 0
                        plot_circle(N(Load(f,1),:), E(s,:), alpha*0.5, dim, load_color, lw);
                    end
                end
            end
    end
end

end

function plot_reaction(N, R, dim, alpha)
assert(any(R(:,1)) <= size(N,1));
r_color = [0.75, 0.75, 0];
lw = 2;
n_Fix = size(R, 1);
circle_r = 0.5;

for f=1:1:n_Fix
    switch dim
        case 2
            quiver(N(R(f,1),1)-alpha*R(f,2),N(R(f,1),2),...
                alpha*R(f,2),0,...
                'Color', r_color, 'LineWidth', lw);
            quiver(N(R(f,1),1),N(R(f,1),2)-alpha*R(f,3),...
                0,alpha*R(f,3),...
                'Color', r_color, 'LineWidth', lw);
            
            if size(R,2)-1 == 3 && R(f,4) ~= 0
                plot_circle(N(R(f,1),:), [], alpha*circle_r, dim, r_color, lw);
            end
        case 3
            quiver3(N(R(f,1),1)-alpha*R(f,2),N(R(f,1),2),N(R(f,1),3),...
                alpha*R(f,2),0,0,...
                'Color', r_color, 'LineWidth', lw);
            quiver3(N(R(f,1),1),N(R(f,1),2)-alpha*R(f,3),N(R(f,1),3),...
                0,alpha*R(f,3),0,...
                'Color', r_color, 'LineWidth', lw);
            quiver3(N(R(f,1),1),N(R(f,1),2),N(R(f,1),3)-alpha*R(f,4),...
                0,0,alpha*R(f,4),...
                'Color', r_color, 'LineWidth', lw);
            
            if size(R,2)-1 == 6
                E = eye(3);
                for s=1:1:3
                    if R(f,4+s) ~= 0
                        plot_circle(N(R(f,1),:), E(s,:), ...
                            alpha, dim, r_color, lw);
                    end
                end
            end
    end
end

end