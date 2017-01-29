% -------------------------------------------------------------------------
% Function:
% label_rm_0 = assign_noise(data, label)
% -------------------------------------------------------------------------
% Target:
% Assign noise sample to the closest cluster
% -------------------------------------------------------------------------
% Input:
% data - m by n data matrix; m: number of samples, n: number of features
% label - m by 1 vector; 
% -------------------------------------------------------------------------
% Output:
% label_rm_0 - m by 1 vector
% -------------------------------------------------------------------------
% Written by Xinyu Zhao
% AU MRI research center, Auburn University
% Janunary, 2017
% https://github.com/xinyuzhao/identification-of-brain-based-disorders

function label_rm_0 = assign_noise(data, label)

k = max(label);
N = size(data, 1);
features = size(data, 2);
M = mean(data);

m = zeros(k, features);
for i = 1 : k
    clusters{i} = data(find(label == i), :);
    m(i, :) = mean(clusters{i}); % Use average as center of cluster
end

label_rm_0 = label;
for i = 1 : N
    distance = [];
    if label(i) == 0
        for j = 1 : k
        distance(j) = sqrt(sum((m(j, :) - data(i, :)) .^ 2));
        [V, I] = min(distance);
        label_rm_0(i) = I;
    end
end

end