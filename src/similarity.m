% -------------------------------------------------------------------------
% Function: 
% sim=similarity(label_1, label_2)
% -------------------------------------------------------------------------
% Target: 
% Compute the cluster similarity between two different clustersing results
% using Torres' method.
% -------------------------------------------------------------------------
% Input: 
% label_1 - m by 1 vector; m: number of samples
% label_2 - m by 1 vector;
% -------------------------------------------------------------------------
% Output: 
% sim - float number; similarity score
% -------------------------------------------------------------------------
% References:
% [1] Torres, G., Basnet, R., Sung, A., 2008. A similarity measure for
% clustering and its applicaitons. Proc. World Acad. Sci. Eng. Technol. 31,
% 490-496.
% -------------------------------------------------------------------------
% Written by Xinyu Zhao
% AU MRI research center, Auburn University
% Janunary, 2017
% https://github.com/xinyuzhao/identification-of-brain-based-disorders

function sim = similarity(label_1, label_2)

cluster_1 = seperate_into_groups(label_1);
num_cluster_1 = length(cluster_1);

cluster_2 = seperate_into_groups(label_2);
num_cluster_2 = length(cluster_2);

C_sim = zeros(num_cluster_1, num_cluster_2);
for i = 1 : num_cluster_1
    for j = 1 : num_cluster_2
        a = cluster_1{i};
        b = cluster_2{j};
        T1 = intersect(a, b);
        T2 = union(a, b);
        C_sim(i, j) = length(T1) / length(T2);
    end
end
sim = sum(sum(C_sim)) / max(num_cluster_1, num_cluster_2);
end




