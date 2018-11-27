clc
clf
close all
format shortE

addpath('utility', 'analysis')

% D2_file_name = 'sf-test_2D_truss.json';
D2_file_name = 'sf-test_2D_beam.json';
% D3_file_name = 'sf-test_3-frame.json';
% D3_file_name = 'cant_2_3D_frame.json';
D3_file_name = 'cant_3D_beam.json';

ins_pth = fullfile(pwd, strcat('test\problem_instances\',D3_file_name));

% input unit:
% nodal coordinates: meter
% pressure (E,G): kN/cm^2
% cross section: cm
[N, T, S, A, m_p] = parse_frame_json(ins_pth);

% Define loads
% L = [node,Qx,Qy,Qz,Mx,My,Mz; ...]
% 2D truss load case
% Load = [3,0,-0.1,0];

% 2D beam
% Load = [2,0,-0.1,0];

Load = [2,0,0,-0.1,0,0,0]; %kN

magnif = 10;

% Output unit: force: kN, length: meter
[F, R, D] = displacement_method(N, T, S, A, m_p, Load, 'Method', 'frame')

% draw_frame(N, T, S, Load, F, R, D, 1, 5, magnif, 0.5);
