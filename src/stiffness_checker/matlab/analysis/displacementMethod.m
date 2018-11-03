% 1.571 - Fall 2018 - Homework 4
% Modified by Yijiang Huang (yijiangh@mit.edu) on 10/10/2018

function [F, R, D] = displacementMethod(N, T, S, L, A, E)

% INPUT:
%
% N = node coordinates
% (number of nodes)-by-2 matrix with N(n,1) and N(n,2) the X and
% Y coordinates of node n
%
% T = truss topology
% (number of elements)-by-2 matrix with T(e,1) and T(e,2) the indices of
% the starting and ending nodes of element e
%
% S = support definition
% (number of fixities)-by-2 matrix with S(s,1) the index of a node and
% S(s,2) = 1 or 2 depending on whether the fixed DOF is in the X or Y
% direction
%
% L = load definition
% (number of loaded nodes)-by-3 matrix with L(n,1) the index of a loaded
% node and L(n,2) and L(n,3) the loads applied to that node in the X and
% Y directions
%
% A = element areas
% (number of elements)-by-1 matrix with A(e,1) the area of element e
%
% E = Young modulus
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

% Program starts here

% Ddetermine the number of nodes and elements
nNodes = size(N,1);
nElements = size(T,1);

% FILL IN: initialize the stiffness matrix K of the structure, with a row
% and a column for each degree of freedom (DOF)
dof = 2 * nNodes;
K = zeros(dof, dof);

id_map = zeros(4,nElements);
for e=1:1:nElements
    id_map(1,e) = T(e,1)*2-1;
    id_map(2,e) = T(e,1)*2;
    
    id_map(3,e) = T(e,2)*2-1;
    id_map(4,e) = T(e,2)*2;
end

for e=1:1:nElements 
    % FILL IN: calculate the terms of the stiffness matrix of element e.
    % Note that terms are repeated in the stiffness matrix of the element,
    % reducing the number of calculations required.
    end_u = N(T(e,1), :);
    end_v = N(T(e,2), :);
    
    Le = norm(end_v - end_u);
    cos = (end_v(1)-end_u(1))/Le;
    sin = (end_v(2)-end_u(2))/Le;
    
    Rot = zeros(2,4);
    Rot(1,1)=cos;
    Rot(1,2)=sin;
    Rot(2,3)=cos;
    Rot(2,4)=sin;
    
    K_loc = (E*A(e,1)/Le) * [[1,-1];[-1,1]];
    Ke = (Rot') * K_loc * Rot;
    
    % FILL IN: add the terms of the stiffness matrix of element e to the
    % stiffness matrix K of the structure
    for i=1:1:4
        l = id_map(i,e);
        for j=1:1:4
            ll = id_map(j,e);
            K(l,ll) = K(l,ll) + Ke(i,j);
        end
    end
    
    % Syntax example: the following command adds the terms of matrix N,
    % assumed 2-by-2, to the terms between rows 3 and 4 and between columns
    % 7 and 8 of a larger matrix M:
    % M(3:4, 7:8) = M(3:4, 7:8) + N;
end

% Assemble the load vector, a matrix of size (number of DOFs)-by-1
Q = zeros(dof,1);
nLoadedNodes = size(L,1);
for i=1:1:nLoadedNodes
    n = L(i,1);
    Qx = L(i,2);
    Qy = L(i,3);
    Q(2*n-1,1) = Qx;
    Q(2*n,1) = Qy;
end

% Assemble the list of fixed DOFs fixedList, a matrix of size
% 1-by-(number of fixities)
nFixities = size(S,1);
fixedList = zeros(1, nFixities);
for f=1:1:nFixities
    n = S(f,1);
    if S(f,2) == 1
        fixedList(f) = 2*n-1;
    else
        fixedList(f) = 2*n;
    end
end

full_f = zeros(1,dof);
for f=1:1:nFixities
    full_f(fixedList(f))=1;
end

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

% FILL IN: create the matrices Kmm and Kfm and the vectors Qm and Qf
% defined in the lecture notes. These matrices and vectors can be
% extracted from K and Q or obtained by copying K and Q and deleting some
% rows and columns from the copies.

Perm = eye(dof);
Perm = Perm(perm_RO,:);

% K * (10.0/(E*A(1,1)))
K_perm = Perm * K * inv(Perm);

Kmm = K_perm(1:nFree,1:nFree);
Kfm = K_perm(nFree+1:dof,1:nFree);

% FILL IN: calculate the vectors Um and Rf defined in the lecture notes.
Q_perm = Perm*Q;
Qm = Q_perm(1:nFree);
Qf = Q_perm(nFree+1:dof);

Um = Kmm\Qm;
Rf = Kfm*Um - Qf;

Us = zeros(nFixities,1);
U_perm = [Um;Us];
U = Perm\U_perm;

Rf_full = [zeros(nFree,1); Rf];
Rf_full = Perm\Rf_full;

for i=1:1:nFixities
    Rf(i) = Rf_full(fixedList(i));
end
% Store the support reactions in the output vector R
R = Rf;

% Reorganize the DOF displacements to form the output matrix D (see output
% description above)
D = zeros(nNodes, 2);
for i=1:1:nNodes
    D(i,1) = U(2*i-1);
    D(i,2) = U(2*i);
end

% Initialize the output vector F, which must be filled with the
% element forces
F = zeros(nElements,1);

for e=1:1:nElements
    % FILL IN: determine the axial force in element e from the
    % displacements of its nodes, and store the results in the output
    % vector F
    end_u = N(T(e,1), :);
    end_v = N(T(e,2), :);
    
    Le = norm(end_v - end_u);
    cos = (end_v(1)-end_u(1))/Le;
    sin = (end_v(2)-end_u(2))/Le;
    
    Rot = zeros(2,4);
    Rot(1,1)=cos;
    Rot(1,2)=sin;
    Rot(2,3)=cos;
    Rot(2,4)=sin;
    
    K_loc = (E*A(e,1)/Le) * [1,-1;-1,1];
    Ue = [U(id_map(1,e));
        U(id_map(2,e));
        U(id_map(3,e));
        U(id_map(4,e))];
    
    Fe = K_loc * Rot * Ue;
    F(e,1) = Fe(2);
end

end