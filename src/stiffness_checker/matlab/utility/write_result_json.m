function write_result_json(abs_file_path, N, T, S, e_F, r_F, n_D)
% all the indices here are converted back to C#/C++'s 0-convention

n_D_data = [];
for i=1:1:size(N,1)
    n_d.node_id = i-1;
    n_d.node_pos = N(i,:);
    n_d.displacement = n_D(i,:);
    n_D_data = [n_D_data, n_d];
end
data.node_displacement = n_D_data;

e_F_data = [];
for i=1:1:size(T,1)
    e_f.element_id = i-1;
    e_f.node_u_id = T(i,1)-1;
    e_f.node_v_id = T(i,2)-1;
    e_f.reaction = e_F(i,:);
    e_F_data = [e_F_data, e_f];
end
data.element_reaction = e_F_data;

r_F_data = [];
for i=1:1:size(S,1)
    r_F_e.node_id = S(i,1)-1;
    r_F_e.node_pos = N(S(i,1),:);
    r_F_e.reaction = r_F(i,:);    
    r_F_data = [r_F_data, r_F_e];
end
data.fixity_reaction = r_F_data;

addpath(fullfile(pwd, 'external\jsonlab-1.5'));
savejson('', data, abs_file_path);

end