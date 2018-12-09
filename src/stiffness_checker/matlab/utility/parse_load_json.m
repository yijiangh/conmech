function [Load, include_self_weight] = parse_load_json(abs_file_path)

addpath(fullfile(pwd, 'external',filesep,'jsonlab-1.5'));
assert(exist(abs_file_path, 'file') == 2);

data = loadjson(abs_file_path);

dim = data.dimension;
n_Loaded_Nodes = length(data.point_load_list);

if dim==2
    node_dof = 3;
else
    node_dof = 6;
end

Load = zeros(n_Loaded_Nodes, 1+node_dof);

for i=1:1:n_Loaded_Nodes
    % C# index starts with 0
    Load(i,1) = data.point_load_list{i}.applied_node_id + 1;
    
    if 2 == dim
        Load(i,2) = data.point_load_list{i}.Fx;
        Load(i,3) = data.point_load_list{i}.Fy;
        Load(i,4) = data.point_load_list{i}.Mz;
    else
        Load(i,2) = data.point_load_list{i}.Fx;
        Load(i,3) = data.point_load_list{i}.Fy;
        Load(i,4) = data.point_load_list{i}.Fz;
        Load(i,5) = data.point_load_list{i}.Mx;
        Load(i,6) = data.point_load_list{i}.My;
        Load(i,7) = data.point_load_list{i}.Mz;
    end
end

include_self_weight = data.include_self_weight;

assert(include_self_weight == 1 || n_Loaded_Nodes > 0);
end