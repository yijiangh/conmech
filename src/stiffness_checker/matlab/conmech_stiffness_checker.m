clc
clf

addpath('utility', 'analysis')

file_name = 'sf-test_3-frame.json';
ins_pth = fullfile(pwd, strcat('test\problem_instances\',file_name));

[N, T, S, A, m_p] = parse_frame_json(ins_pth);

% Define loads
% L = [node,Qx,Qy,Qz,Mx,My,Mz; ...]
L = [4,0,0,-100,0,0,0]; % (kipf)

% [F, R, D] = displacementMethod(N, T, S, L, A, E)
magnif = 10;
draw_frame(N, T, S, [], [], 1, 5, magnif);