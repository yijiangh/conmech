function [N, T, S, A, m_p] = parse_frame_json(abs_file_path)
% Inputs:
%   abs_file_path: absolute file path.
% Outputs:
%   N: #nodes x 3 nodal position list. 
%       [x,y,z; ...]
%   T: #elements x 2 element connectivity list. 
%       [start_node_id, end_node_id; ...]
%   S: #fix x 2 fixities list.
%       [node_id, direction[1:6]; ...],
%       direction = 1-3: translation x,y,z
%       direction = 4-6: rotation Rxx, Ryy, Rzz
%   A: #elements x 1 cross section list.
%       [A_i;...]
%   m_p: a struct containing material properties
%       E: MPa
%       G: MPa
%       gamma: [unit_less]
%       density: kg/m3

addpath(fullfile(pwd, 'external\jsonlab-1.5'));
assert(exist(abs_file_path, 'file') == 2);

data = loadjson(abs_file_path);

n_Nodes = length(data.node_list);
n_Elements = length(data.element_list);

N = zeros(n_Nodes, 3);
T = zeros(n_Elements, 2);

% construct nodal and fixities list
S = [];
for i=1:1:n_Nodes
    N(i,1) = data.node_list{i}.point.X;
    N(i,2) = data.node_list{i}.point.Y;
    N(i,3) = data.node_list{i}.point.Z;
    
    if data.node_list{i}.is_grounded
        S = [S; [i, linspace(1,1,6)]];
    end
end

% construct element list
for i=1:1:n_Elements
    % initial file's nodal id starts from 0
    T(i,:) = data.element_list{i}.end_node_ids + 1;
end

assert(0 == any(any(T > n_Nodes) ~= 0));

m_p.E = data.material_properties.youngs_modulus;
m_p.G = data.material_properties.shear_modulus;
m_p.gamma = data.material_properties.poisson_ratio;
m_p.density = data.material_properties.density;

% construct cross sections
r = data.material_properties.radius;
A = ones(n_Elements,1) * pi * r^2;

end