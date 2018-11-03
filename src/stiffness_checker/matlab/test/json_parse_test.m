% test file for json file parsing using jsonlab
% https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab-a-toolbox-to-encode-decode-json-files

addpath(fullfile(pwd, 'external\jsonlab-1.5'));
test_file = fullfile(pwd, 'test\problem_instances\sf-test_3-frame.json');
assert(exist(test_file, 'file') == 2);

data = loadjson(test_file);
data.element_list{1}.end_node_ids