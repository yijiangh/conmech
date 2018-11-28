function [F, RF, D] = displacement_method(N, T, S, m_p, Load, self_weight, varargin)
% INPUT:
%
% N = node coordinates
% (number of nodes)-by-dim matrix, [x,y;x,y; ...] or [x,y,z; ...]
%
% T = truss topology
% (number of elements)-by-2 matrix with T(e,1) and T(e,2) the indices of
% the starting and ending nodes of element e
%
% S = support definition
% (number of fixities)-by-(4 or 7) matrix with S(s,1) the index of a node and
% S(s,2:4 or 7) = [Tx, Ty, Tz, Rx, Ry, Rz] (boolean)
%
% A = element areas
% (number of elements)-by-1 matrix with A(e,1) the area of element e
%
% m_p: a struct containing various material properties
%
% Load = load definition
% (number of loaded nodes)-by-(4 or 7) matrix with L(n,1) the index of a loaded
% node and L(n,2) and L(n,3) the loads applied to that node in the X and
% Y directions
%
% dim = dimension (2 or 3)
%
% OUTPUT:
%
% F = element forces
% (number of elements)-by-1 matrix with F(e,1) the force in element e
%
% R = support reactions
% (number of fixities)-by-1 matrix with R(f,1) the reaction developed by
% fixity f
%
% D = node displacements
% (number of nodes)-by-2 matrix with D(n,1) and D(n,2) the displacements of
% node n in the X and Y directions
%
% Note: the calculation inside use centimeter and kN
% output will be in meter and kN

% Ddetermine the number of nodes and elements
nNodes = size(N,1);
dim = size(N,2);
nElements = size(T,1);

% check dimension of S
if dim == 2
    assert(4 == size(S,2));
else
    assert(7 == size(S,2));
end

% default values
% if 2 == dim
%     method = 'truss';
% else
%     method = 'frame';
% end

% Map of parameter names to variable names
params_to_variables = containers.Map({'Method'},{'method'});
a_m = available_methods();
v = 1;
while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
        assert(v+1<=numel(varargin));
        v = v+1;
        % Trick: use feval on anonymous function to use assignin to this workspace
        % see: http://www.alecjacobson.com/weblog/?p=3792
        feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
        error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
end
if isempty(find(strcmp(a_m,method), 1))
    error('Unsupported displacement method: %s', method);
end

% initialize the stiffness matrix K of the structure, with a row
% and a column for each degree of freedom (DOF)
if dim == 2
    switch method
        case 'truss'
            node_dof = 2;
        case 'frame'
            node_dof = 3;
    end
    full_node_dof = 3;
else
    switch method
        case 'truss'
            node_dof = 3;
        case 'frame'
            node_dof = 6;
    end
    full_node_dof = 6;
end

dof = node_dof * nNodes;
K = zeros(dof, dof);

id_map = zeros(nElements,node_dof*2);
dof_lin = linspace(node_dof-1,0,node_dof);
for e=1:1:nElements
    % end nodal id
    u = T(e,1);
    v = T(e,2);
    
    id_map(e, 1:node_dof) = u*node_dof*linspace(1,1,node_dof)-dof_lin;
    id_map(e, node_dof+1:2*node_dof) = v*node_dof*linspace(1,1,node_dof)-dof_lin;
end

% assuming uniform cross sections now
[Jx, Iy, Iz] = cross_sec_properties(m_p);
A = ones(nElements,1) * pi * m_p.r^2;

% roll anlgle for local frame
roll_angle = 0;

switch method
    case 'truss'
        if 2 == dim
            ex_id = [1,4];
            xy_id = [1,2,4,5];
        end
        if 3 == dim
            ex_id = [1,7];
            xy_id = [1,2,3,7,8,9];
        end
    case 'frame'
        ex_id = linspace(1,node_dof*2,node_dof*2);
        xy_id = linspace(1,node_dof*2,node_dof*2);
end
e_react_dof = size(ex_id,2)/2;

K_loc_list = {};
R_list = {};
for e=1:1:nElements
    % calculate the terms of the stiffness matrix of element e.
    % Note that terms are repeated in the stiffness matrix of the element,
    % reducing the number of calculations required.
    
    end_u = N(T(e,1), :);
    end_v = N(T(e,2), :);
    
    R_b = local_frame(end_u, end_v, roll_angle);
    K_loc = local_stiffness_matrix(norm(end_u-end_v), A(e),...
        Jx, Iy, Iz, dim, m_p);
    
    R = zeros(full_node_dof*2, full_node_dof*2);
    for k=1:1:(full_node_dof/3)*2
        R(k*3-3+1:k*3,...
          k*3-3+1:k*3) = R_b;
    end
    
    R = R(ex_id, xy_id);
    K_loc = K_loc(ex_id, ex_id);
    
    K_loc_list{end+1} = K_loc;
    R_list{end+1} = R;
    
    Ke = (R') * K_loc * R;
    
    % add the terms of the stiffness matrix of element e to the
    % stiffness matrix K of the structure
    for i=1:1:(2*node_dof)
        s = id_map(e,i);
        for j=1:1:(2*node_dof)
            t = id_map(e,j);
            K(s,t) = K(s,t) + Ke(i,j);
        end
    end
end

% Assemble the load vector, a matrix of size (number of DOFs)-by-1
Q = zeros(dof,1);
nLoadedNodes = size(Load,1);
if ~isempty(Load)
    assert(size(Load,2)==1+node_dof);
    for i=1:1:nLoadedNodes
        n = Load(i,1);
        q_id = linspace(1,node_dof,node_dof);
        Q_node = Load(i, 1 + q_id);
        Q(node_dof*n-node_dof+1 : node_dof*n,1) = Q_node;
    end
end

if self_weight
    Q_sw = zeros(dof,1);
    
    for e=1:1:nElements
        % density is in kN/m3, radius in cm, length in m
        end_u = N(T(e,1), :);
        end_v = N(T(e,2), :);
        Le = norm(end_v-end_u);
        q_sw = A(e)*m_p.density;
                
        if method == 'frame'
            if dim == 2
                R_eG = R_list{e};
                
                Q_sw_G = zeros(2*node_dof,1);
        
                % end gravity in global frame
                Q_sw_G(dim) = -q_sw*Le/2;
                Q_sw_G(node_dof+dim) = -q_sw*Le/2; 

                % force density in local frame
                Q_sw_e = R_eG*Q_sw_G;
                fixed_end_sw_load = Q_sw_e;
                fixed_end_sw_load(3) = Q_sw_e(2)/Le * Le^2 / 12;
                fixed_end_sw_load(6) = -fixed_end_sw_load(3);
                fixed_end_sw_load = R_eG' * fixed_end_sw_load;
            else
%                 fixed_end_sw_load(5) = Q_sw_e(2)/Le * Le^2 / 12;
%                 fixed_end_sw_load(6) = Q_sw_e(3)/Le * Le^2 / 12;
%                 fixed_end_sw_load(11) = -fixed_end_sw_load(5);
%                 fixed_end_sw_load(12) = -fixed_end_sw_load(6);
%                 fixed_end_sw_load = R_eG' * fixed_end_sw_load;

                % https://github.mit.edu/yijiangh/frame3dd/blob/master/frame3dd_io.c#L905
                R_eG = local_frame(end_u, end_v, roll_angle);
                fixed_end_sw_load = zeros(2*node_dof,1);
                
                fixed_end_sw_load(3) = -q_sw * Le / 2.0;
                fixed_end_sw_load(4) = q_sw * Le^2 / 12.0 *...
                    ((-R_eG(2,1)*R_eG(3,3)+R_eG(2,3)*R_eG(3,1))*(-1));
                fixed_end_sw_load(5) = q_sw * Le^2 / 12.0 *...
                    ((-R_eG(2,2)*R_eG(3,3)+R_eG(2,3)*R_eG(3,2))*(-1));
                fixed_end_sw_load(6) = 0;
                
                fixed_end_sw_load(9) = -q_sw * Le / 2.0;
                fixed_end_sw_load(10) = -q_sw * Le^2 / 12.0 *...
                    ((-R_eG(2,1)*R_eG(3,3)+R_eG(2,3)*R_eG(3,1))*(-1));
                fixed_end_sw_load(11) = -q_sw * Le^2 / 12.0 *...
                    ((-R_eG(2,2)*R_eG(3,3)+R_eG(2,3)*R_eG(3,2))*(-1));
                fixed_end_sw_load(12) = 0;
            end
        end
        
        Q_sw(id_map(e,1:end)) = Q_sw(id_map(e,1:end)) + fixed_end_sw_load;
    end
    
    Q = Q + Q_sw;
end

% Assemble the list of fixed DOFs fixedList, a matrix of size
% 1-by-(number of fixities)
S = S(:,1:1+node_dof);
nFixities = sum(sum(S(:,2:1+node_dof)==1));
full_f = zeros(1,dof);
for f=1:1:size(S,1)
    n = S(f,1);
    full_f(n*node_dof-node_dof+1 : n*node_dof)=S(f,2:1+node_dof);
end
assert(sum(full_f(:)==1) == nFixities);
% fixedList = find(full_f==1);

nFree = dof - nFixities;
free_tail=1;
fix_tail=nFree+1;
% generate perm id map
perm_RO = 1:1:dof;
for i=1:1:dof
    if 0 == full_f(i)
        perm_RO(free_tail)=i;
        free_tail = free_tail + 1;
    else
        perm_RO(fix_tail)=i;
        fix_tail = fix_tail + 1;
    end
end
assert(size(perm_RO,2) == size(unique(perm_RO),2));

% create the matrices Kmm and Kfm and the vectors Qm and Qf
% defined in the lecture notes. These matrices and vectors can be
% extracted from K and Q or obtained by copying K and Q and deleting some
% rows and columns from the copies.

Perm = eye(dof);
Perm = Perm(perm_RO,:);

K_perm = Perm * K * inv(Perm);

Kmm = K_perm(1:nFree,1:nFree);
Kfm = K_perm(nFree+1:dof,1:nFree);

% FILL IN: calculate the vectors Um and Rf defined in the lecture notes.
Q_perm = Perm*Q;
Qm = Q_perm(1:nFree);
Qf = Q_perm(nFree+1:dof);

Um = Kmm\Qm; % in m
Rf = Kfm*Um - Qf;

Us = zeros(nFixities,1);
U_perm = [Um;Us];
U = Perm\U_perm;

Rf_full = [zeros(nFree,1); Rf];
Rf_full = Perm\Rf_full;

RF = S(:,2:end);
for f=1:1:size(S,1)
    fix_node_id = S(f,1);
    RF(f,:) = Rf_full(fix_node_id*node_dof-node_dof+1 : fix_node_id*node_dof);
end

% Reorganize the DOF displacements to form the output matrix D  (see output
% description above)
D = zeros(nNodes, node_dof);
for i=1:1:nNodes
    D(i,1:node_dof)=U(i*node_dof-node_dof+1:i*node_dof);
end

% Initialize the output vector F, which must be filled with the
% element forces
F = zeros(nElements,e_react_dof*2);

for e=1:1:nElements
    % determine the axial force in element e from the
    % displacements of its nodes, and store the results in the output
    % vector F
    Ue = U(id_map(e,1:end));
    
    Fe = K_loc_list{e} * R_list{e} * Ue;
    F(e,:) = Fe';
end

end
