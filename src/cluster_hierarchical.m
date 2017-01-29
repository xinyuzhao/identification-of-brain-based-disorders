% -------------------------------------------------------------------------
% Function:
% [T, V] = cluster_hierarchical(data, max_num)
% -------------------------------------------------------------------------
% Target:
% Cluster data using hierarchical clustering
% -------------------------------------------------------------------------
% Input:
% data - m by n data matrix; m: number of samples, n: number of
% max_num - integer; maximum number of clusters
% -------------------------------------------------------------------------
% Output:
% T - m by 1 vector
% V - float number; clustering score
% -------------------------------------------------------------------------
% References: 
% [1] Liao, W., Chen, H., Yang, Q., Lei, X., 2008. Analysis of fMRI data
% using improved self-organizing mapping and spatio-temporal metric
% hierarchical clustering. IEEE Trans. Med. Imaging 27, 1472-1483.
% -------------------------------------------------------------------------
% Written by Xinyu Zhao
% AU MRI research center, Auburn University
% Janunary, 2017
% https://github.com/xinyuzhao/identification-of-brain-based-disorders

function [T, V] = cluster_hierarchical(data, max_num)

distM = pdist(data);
Z = linkage(distM, 'average');

CH_indice = zeros(max_num, 1);
for num_c = 2 : max_num
    T = cluster(Z, 'maxclust', num_c);
    CH_indice(num_c) = compute_CH_index(data, T);
end
[V, I] = max(CH_indice);
T = cluster(Z, 'maxclust', I);
end


