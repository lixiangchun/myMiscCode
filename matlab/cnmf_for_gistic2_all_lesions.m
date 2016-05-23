

function cnmf_for_gistic2_all_lesions(all_lesions_file)

%
%
% all_lessions_file: gistic2 focal CNA file
%
%

a = dataset('File', all_lesions_file);
n = length(a.AmplitudeThreshold);
if mod(n, 2) ~= 0
	error('mod(nrow_of_input_file, 2) ~= 0, ')
end

[nr, nc] = size(a);
a = a((n/2 + 1):nr, 10:(nc - 1));

a = double(a);

kstart = 2;
kend = 11; 
nloop = 50;
verbose = 1;
[pathstr, prefix, ext] = fileparts(all_lesions_file);
cnmf(a, kstart, kend, nloop, verbose, prefix);

