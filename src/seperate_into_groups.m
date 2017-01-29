% -------------------------------------------------------------------------
% Function: 
% clusters=seperate_into_groups(label)
% -------------------------------------------------------------------------
% Target: convert a vector into cell array
% -------------------------------------------------------------------------
% Input: 
% label - m by 1 vector; m: number of samples
% -------------------------------------------------------------------------
% Output: 
% clusters - cell array; the size of cell array equals to the number of
% unique numbers in label
% -------------------------------------------------------------------------
% Written by Xinyu Zhao
% AU MRI research center, Auburn University
% Janunary, 2017
% https://github.com/xinyuzhao/identification-of-brain-based-disorders

function clusters = seperate_into_groups(label)
x_label = label;

row = size(x_label, 1);
data_index = 1 : row;
I = unique(x_label);

clusters = [];
for i = 1 : length(I)
    clusters{i} = data_index(1, find(x_label == I(i)));
end

clusters = clusters(~cellfun('isempty', clusters));
end




