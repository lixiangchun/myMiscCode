

function cnmf(a, kstart, kend, nloop, verbose, prefix)
%
% Consensus clustering based on nonnegative matrix factorization (NMF)
%
% Usage: cnmf(a, kstart, kend, nloop, verbose, prefix)
%
% Parameters:
%       a:       input matrix, for example, genes-by-samples
%       kstart:  start NMF rank
%       kend:    end NMF rank
%       verbose: 1 - verbose, 0 - silence
%       prefix:  prefix of output files
%
% This function produces 4 output files:
% `prefix`.mat:  saved ordered consensus matrix, cluster ID, cophenetic and dispersion coeffs.
% `prefix`.clustid.txt: output file of cluster ID.
% `prefix`.pdf:  consensus clustering heatmap in pdf format.
% `prefix`.fig:  consensus clustering heatmap in fig format.
%
% The NMF algorithm used is nnmf in Matlab.
%
% Reference: Brunet et al. Metagenes and molecular pattern discovery using matrix factorization.
%

tic

if sum(any(a < 0)) > 0
	a=[max(a,0);-min(a,0)];
end

consensus = nmfconsensus(a,kstart,kend,nloop,verbose);
[ordcons,clustid,ordindex,coph] = nmforderconsensus(consensus,kstart,kend);
rho = dispersionCoefficients(consensus, kstart, kend);

save([prefix,'.mat'], 'ordcons', 'clustid', 'coph', 'rho');
dlmwrite([prefix, '.clustid.txt'], clustid, '\t');

figure; clf;
plotnmforderconsensus(ordcons, coph, rho)
fig = gcf;
%fig.PaperSize = [8, 5]; % works for matlab2015b
%fig.PaperPosition = [0.15, 0.4, 8, 4]; % works for matlab2015b
set(fig, 'PaperSize', [8, 5]);
set(fig, 'PaperPosition', [0.15, 0.4, 8, 4]);

%% default colormap in matlab2015b
parula = [0.2081 0.1663 0.5292
0.2116 0.1898 0.5777
0.2123 0.2138 0.6270
0.2081 0.2386 0.6771
0.1959 0.2645 0.7279
0.1707 0.2919 0.7792
0.1253 0.3242 0.8303
0.0591 0.3598 0.8683
0.0117 0.3875 0.8820
0.0060 0.4086 0.8828
0.0165 0.4266 0.8786
0.0329 0.4430 0.8720
0.0498 0.4586 0.8641
0.0629 0.4737 0.8554
0.0723 0.4887 0.8467
0.0779 0.5040 0.8384
0.0793 0.5200 0.8312
0.0749 0.5375 0.8263
0.0641 0.5570 0.8240
0.0488 0.5772 0.8228
0.0343 0.5966 0.8199
0.0265 0.6137 0.8135
0.0239 0.6287 0.8038
0.0231 0.6418 0.7913
0.0228 0.6535 0.7768
0.0267 0.6642 0.7607
0.0384 0.6743 0.7436
0.0590 0.6838 0.7254
0.0843 0.6928 0.7062
0.1133 0.7015 0.6859
0.1453 0.7098 0.6646
0.1801 0.7177 0.6424
0.2178 0.7250 0.6193
0.2586 0.7317 0.5954
0.3022 0.7376 0.5712
0.3482 0.7424 0.5473
0.3953 0.7459 0.5244
0.4420 0.7481 0.5033
0.4871 0.7491 0.4840
0.5300 0.7491 0.4661
0.5709 0.7485 0.4494
0.6099 0.7473 0.4337
0.6473 0.7456 0.4188
0.6834 0.7435 0.4044
0.7184 0.7411 0.3905
0.7525 0.7384 0.3768
0.7858 0.7356 0.3633
0.8185 0.7327 0.3498
0.8507 0.7299 0.3360
0.8824 0.7274 0.3217
0.9139 0.7258 0.3063
0.9450 0.7261 0.2886
0.9739 0.7314 0.2666
0.9938 0.7455 0.2403
0.9990 0.7653 0.2164
0.9955 0.7861 0.1967
0.9880 0.8066 0.1794
0.9789 0.8271 0.1633
0.9697 0.8481 0.1475
0.9626 0.8705 0.1309
0.9589 0.8949 0.1132
0.9598 0.9218 0.0948
0.9661 0.9514 0.0755
0.9763 0.9831 0.0538];
%fig.Colormap = parula; % works for matlab2015b
set(fig, 'Colormap', parula);
saveas(fig, [prefix,'.pdf'], 'pdf');
saveas(fig, [prefix,'.fig'], 'fig');
close all;

toc


function  [ordcons,clustid,ordindex,coph] = nmforderconsensus(consensus,kstart,kend)
%
% Jean-Philippe Brunet
% Cancer Genomics
% The Broad Institute
% brunet@broad.mit.edu
%
% This software and its documentation are copyright 2004 by the
% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
% or functionality. 
%
%
% Reordering of consensus matrices using hierarchical clustering
%
% consensus : 3d array of consensus matrices 
%             dimensions : kend x M x M 
%
% kstart, kend : range of  consensus matrices to reorder
%
%
% ordcons  : 3d array of consensus matrices 
%            dimensions : kend x M x M
%            Values in ordcons(1:kstart-1,:,:) should be ignored
%
% clustid : dimensiosn kend x M
%           clustid(k,:) has cluster id for each sample (in original order)
%           Values in clustid(1:kstart-1,:) should be ignored
%
% ordindex : dimension kend x M 
%            ordindex (k,:) is the permutation vector 
%            for the kth consensus matrix 
%            Values in ordindex(1:kstart-1,:) should be ignored
%
% coph :  dimension kend
%         coph(k) has cophenetic coefficients for k 
%         Values in coph(1:kstart-1) should be ignored
%


[kmax,m,m]=size(consensus);

%check that kend is not out of range
if(kend>kmax) 
error('kend is out of range');
end

ordcons=zeros(kmax,m,m);
ordindex=zeros(m,kmax);
clustid=zeros(m,kmax);
coph=zeros(1,kmax);

for i=kstart:kend
u=reshape(consensus(i,:,:),m,m);
[ordcons(i,:,:),clustid(:,i),ordindex(:,i),coph(i)] = nmforderconsensus0(u,i);
end



function  [aordered,clust,ord,coph] = nmforderconsensus0(a,k)
%
% Jean-Philippe Brunet
% Cancer Genomics 
% 

[n,m]=size(a);
ordl=zeros(1,m);
incr=1;

uvec=a(1,2:end);

for i=2:n-1;
uvec=[uvec a(i,i+1:end)]; %get upper diagonal elements of consensus
end

y=1-uvec;                 % consensus are similarities, convert to distances
z=linkage(y,'average');   % use average linkage
coph=cophenet(z,y);

fig = figure('visible','off'); % turn off dendrogram plot
[h,t,ord]=dendrogram(z,0); % get permutation vector
close(fig)

clust=cluster(z,k);       % get cluster id 
aordered=a(ord,ord);
ord=ord';




function consensus = nmfconsensus(a,kstart,kend,nloop,verbose)
%
% Jean-Philippe Brunet
% Cancer Genomics 
% The Broad Institute
% brunet@broad.mit.edu
%
% This software and its documentation are copyright 2004 by the
% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
% or functionality. 
%
% Model selection for NMF
%
% a (n,m) : N (genes) x M (samples) original matrix 
%           numerical data only. Must be positive.
% 
% kstart, kend : range of values of k to test consensus for.
%
% nloop : number of initial conditions per k 
%         (start with 10 or 20, check with more)
%
% verbose : prints iteration count and changes in connectivity matrix elements
%           if not set to 0 
%
% consensus : 3d array of consensus matrices 
%             dimensions : kend x M x M 
%             Values in consensus(1:kstart-1,:,:) should be ignored
%              

% test for negative values in v
if min(min(a)) < 0
error('matrix entries can not be negative');
return
end
if min(sum(a,2)) == 0
error('not all entries in a row can be zero');
return
end

[n,m]=size(a);

consensus=zeros(kend,m,m);
conn=zeros(m,m); 
incr=0;

for j=kstart:kend

if verbose fprintf(1,'rank %d\n',j), end

connac=zeros(m,m); 

for iloop=1:nloop;

incr=incr+1;

if verbose fprintf(1,' iteration %d\n',iloop), end 

%[w,h]=nmf(a,j,verbose);

options = statset('TolX', 1e-8, 'TolFun', 1e-8, 'MaxIter', 100000000);
[w, h] = nnmf(a, j, 'algorithm', 'als', 'replicates', 50, 'options', options); 

%
% compute consensus matrix
%
conn=nmfconnectivity(h); 
connac=connac+conn; % accumulate connectivity matrices

end

consensus(j, :, :)=connac/nloop; %average

end



function conn= nmfconnectivity(h)
%
% Jean-Philippe Brunet
% Cancer Genomics 6/10/03
%
mm=size(h);
k=mm(1);
m=mm(2);


% compute m x m matrix which is 1 if samples are together, 0 elsewhere

% determine sample assignment by its largest metagene expresion value
[y,index]=max(h,[],1); 

mat1=repmat(index,m,1); % spread index down
mat2=repmat(index',1,m); % spread index right

conn=mat1==mat2; % 1 when for pair of samples with same assignement



function plotnmforderconsensus(cons, coph, rho)

k0 = 2;
% 
% plotting cNMF consensus matrices
% Usage: plotnmforderconsensus(ordcons, coph, rho)
%

[k,m,m] = size(cons);
j = 1;
nr = 3;
nc = 4;

if size(cons,1) > nr * nc
    nr = sqrt(size(cons,1) + 3);
    nc = nr;
end

for i = 1:(k-1)
	subplot(nr, nc,i)
	%subplot(size(cons,1)/2,2,i)
	u = cons(k0+i-1,:,:);
	imagesc(squeeze(u));
	title(sprintf('k = %d', k0+i-1));
    j = i;
end

subplot(nr,nc,j+1)
plot(coph, '-o')
xlabel('cNMF rank')
ylabel('Cophenetic')
%ylim([0.7,1])
xlim([2, size(coph,2)])
grid on

subplot(nr,nc,j+2)
plot(rho, '-o')
xlabel('cNMF rank')
%ylabel('\rho')
%ylabel('\tau')
ylabel('Dispersion')
%ylim([0.7,1])
xlim([2, size(rho,2)])
grid on

%disp('Finished calling consplot.')



function [w, h] = normalizewh(w0, h0)
%
% Normalize input matrices w0 and h0 return by nmf(...). The column sums of
% returned w all equal 1, and column sums of h approximately equal column 
% sums of matrix a as input for nmf(...), i.e. sum(w) == 1 and sum(h) == sum(a).
%
% Usage: [w, h] = normalizewh(w0, h0)
%
% Example:
% load iris_dataset
% [w0, h0]=nmf(irisInputs, 2, 0);
% [w, h] = normalizewh(w0, h0);

% source code from mutational signature deciphering framework - extract.m

% initialize w and h
w = rand(size(w0));
h = rand(size(h0));

% the total number of mutational processes
totalProcesses = size(w0, 2);

for j = 1 : totalProcesses
	total = sum( w0(:, j) );
	w(:, j) = w0(:, j) / total;
	h(j, :) = h0(j, :) * total;
end



function rho = dispersionCoefficient(con)

%
% con: consensus matrix for a given rank
%

tot = 0.0;
[n, m] = size(con);
for i = 1:n
	for j = 1:n
		tot = tot + 4.0 * (con(i,j) - 0.5)^2;
	end
end
rho = tot / (n * n);	


function rhos = dispersionCoefficients(cons, kstart, kend)

% Hyunsoo Kim and Haesun Park  Sparse non-negative matrix factorizations 
%+via alternating non-negativity-constrained least squares for microarray 
%+data analysis Bioinformatics 2007(23):1495-1502.

% Input parameter:
% cons: Consensus matrices returned from nmforderconsensus(...)
% kstart, kend: the start and end of NMF rank
%
% Output:
%   rhos      An array of dispersion coefficients.
%             Values in ordindex(1:kstart-1,:) should be ignored
%

rhos = repmat(0,1,kend);
for i = kstart: kend
	con = squeeze(cons(i,:,:));
	rhos(i) = dispersionCoefficient(con);
end
