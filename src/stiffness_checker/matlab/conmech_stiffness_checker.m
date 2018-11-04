clc
clf
close all

addpath('utility', 'analysis')

D2_file_name = '2D_truss_ft.json';
D3_file_name = 'sf-test_3-frame.json';

ins_pth = fullfile(pwd, strcat('test\problem_instances\',D3_file_name));

% input unit:
% nodal coordinates: meter
% pressure (E,G): kN/cm^2
% cross section: cm
[N, T, S, A, m_p] = parse_frame_json(ins_pth);

% Define loads
% L = [node,Qx,Qy,Qz,Mx,My,Mz; ...]
%Load = [3,0,-100,0];
Load = [4,0,0,-0.1,0,0,0]; %kN

magnif = 100;
% draw_frame(N, T, S, [], [], [], 0, 1, 5, magnif);

% Output unit: force: kN, length: meter
[F, R, D] = displacement_method(N, T, S, A, m_p, Load, 'Method', 'truss')

draw_frame(N, T, S, Load, F, R, D, 1, 5, magnif, 5);