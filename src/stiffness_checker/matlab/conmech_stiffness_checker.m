clc
clf
close all
format shortE

addpath('utility', 'analysis')

% frame_file_name = 'sf-test_2D_truss.json';
% frame_file_name = 'sf-test_2D_beam.json';
% frame_file_name = '2D_frame.json';

% frame_file_name = 'sf-test_3-frame.json';
% frame_file_name = 'sf-test_3-frame_0.json';

% frame_file_name = 'cant_2_3D_frame.json';
% frame_file_name = 'cant_3D_beam.json';
% frame_file_name = 'topopt100_3D.json';
frame_file_name = 'tower_3D.json';

ins_pth = fullfile(pwd, strcat('test', filesep, 'problem_instances', filesep, frame_file_name));

% load_file_name = 'topopt100_3D_load_case.json';
load_file_name = 'tower_3D_load_case.json';
load_pth = fullfile(pwd, strcat('test',filesep,'problem_instances',filesep,load_file_name));

% input unit:
% nodal coordinates: meter
% pressure (E,G): kN/m^2
% cross section: m
use_self_weight = 0;
[N, T, S, m_p] = parse_frame_json(ins_pth);
[Load, use_self_weight] = parse_load_json(load_pth);

% Define loads
% L = [node,Qx,Qy,Qz,Mx,My,Mz; ...]

% 2D truss load case
% Load = [3,0,-0.1,0];

% 2D beam
% Load = [2,0,-0.1,0];

% 2D frame
% Load = [4,0,-1,0];

% 3D cases
% 3D beam, node 1 or 2
% Load = [2, 0, 0, -1, 0,0,0]; %kN

% 3D frame
% Load = [1, 0,0,-1, 0,0,0]; %kN
% use_self_weight = 1;

% if(use_self_weight)
%     Load = [];
% end

magnif = 100;

% Output unit: force: kN, length: meter
[element_F, reaction_F, nodal_displ] = displacement_method(N, T, S, m_p, Load, use_self_weight, 'Method', 'frame')

draw_frame(N, T, S, Load, element_F, reaction_F, nodal_displ, 1, 5, magnif, 0.1);

% result_file_name = 'tower_3D_result.json';
% result_save_path = fullfile(pwd, strcat('test',filesep,'cm_results',filesep,result_file_name));
% 
% write_result_json(result_save_path, N, T, S, element_F, reaction_F, nodal_displ)

% % sort displacement
% s_D = zeros(size(nodal_displ,1),1+size(nodal_displ,2));
% for i=1:1:size(s_D,1)
%    S_D(i,1) = norm(nodal_displ(i,:));
% %    nodal_displ(i,:)
%    S_D(i,2:7) = nodal_displ(i,:);
% end
% 
% s_D = sortrows(s_D);
% S_D(end-10:end-5, 2:7)

