

function cnmf_for_broad_values_by_arm(broad_values_by_arm_file)

%
%
% broad_values_by_arm_file: gistic2 broad_values_by_arm.txt
%
%

a = dataset('File', broad_values_by_arm_file);
[nr, nc] = size(a);

a = double(a(:,2:size(a,2)));
input_matrix_dimensions = size(a)

kstart = 2;
kend = 11; 
nloop = 50;
verbose = 1;
[pathstr, prefix, ext] = fileparts(broad_values_by_arm_file);
cnmf(a, kstart, kend, nloop, verbose, prefix);

