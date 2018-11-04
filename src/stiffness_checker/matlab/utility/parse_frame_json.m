function [N, T, S, A, m_p] = parse_frame_json(abs_file_path)
% Inputs:
%   abs_file_path: absolute file path.
% Outputs:
%   N: #nodes x 3 nodal position list. (meter)
%       [x,y,z; ...]
%   T: #elements x 2 element connectivity list. 
%       [start_node_id, end_node_id; ...]
%   S: #fix x 2 fixities list.
%       [node_id, direction[1:6]; ...],
%       direction = 1-3: translation x,y,z
%       direction = 4-6: rotation Rxx, Ryy, Rzz
%   A: #elements x 1 cross section list. (mm^2)
%       [A_i;...]
%   m_p: a struct containing material properties
%       E (Young's Modulus): MPa
%       G (Shear Modulus): MPa
%       gamma (poisson ratio): [unit_less]
%       density: kg/m3

addpath(fullfile(pwd, 'external\jsonlab-1.5'));
assert(exist(abs_file_path, 'file') == 2);

data = loadjson(abs_file_path);

dim = data.dimension;
n_Nodes = length(data.node_list);
n_Elements = length(data.element_list);

N = zeros(n_Nodes, dim);
T = zeros(n_Elements, 2);

% construct nodal and fixities list
S = [];
for i=1:1:n_Nodes
    N(i,1) = data.node_list{i}.point.X;
    N(i,2) = data.node_list{i}.point.Y;
    if 3 == dim
        N(i,3) = data.node_list{i}.point.Z;
    end
    
    if data.node_list{i}.is_grounded
        if 2 == dim
            assert(3 == length(data.node_list{i}.fixities));
        else
            assert(6 == length(data.node_list{i}.fixities));
        end
        
        S = [S; [i, data.node_list{i}.fixities]];
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
m_p.mu = data.material_properties.poisson_ratio;
m_p.density = data.material_properties.density;
m_p.r = data.material_properties.radius;

% construct cross sections
A = ones(n_Elements,1) * pi * m_p.r^2;

end