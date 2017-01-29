% -------------------------------------------------------------------------
% Function:
% [x_label, V] = clustering_optics(x, num_c)
% -------------------------------------------------------------------------
% Target:
% Cluster data using ordering points to identify the clustering structures.
% -------------------------------------------------------------------------
% Input:
% x - m by n data matrix; m: number of samples, n: number of
% num_c - integer; maximum number of clusters
% -------------------------------------------------------------------------
% Output:
% x_label - m by 1 vector
% V - float number; clustering score
% -------------------------------------------------------------------------
% References: 
% [1] Ankerst, M., Breuning, M.M., Kriegel, H.-P., Sander, J., 1999.
% Optics: Ordering points to identify the clustering structures. ACM Sigmod
% Rec., 49-60.
% -------------------------------------------------------------------------
% Written by Xinyu Zhao
% AU MRI research center, Auburn University
% Janunary, 2017
% https://github.com/xinyuzhao/identification-of-brain-based-disorders

function [x_label, V] = cluster_optics(x, num_c)
[row, col] = size(x);
[RD,CD,order]=optics(x, ceil(row / 25)); % MinPts = data_size * 0.5%
Y_range = max(max(RD), max(CD));
X_range = size(x, 1);
data_size = X_range;

[M, N] = size(x);

VRC = [];
threshold = [];
num = 1;
for ei = max(RD) : -mean(RD) / 20 : mean(RD) * 0.5
    clusterId = 0; % 0 represents noise
    x_label = zeros(data_size, 1);
    for i = order
        if RD(1, i) > ei
            if CD(1, i) <= ei
                clusterId = clusterId + 1;
                x_label(i) = clusterId;
            else
                x_label(i) = 0;
            end
        else
            x_label(i) = clusterId;
        end
    end
    
    %x_label = assign_noise(x, x_label');
    
    if max(x_label) <= 1 || max(x_label) > num_c
        VRC(num) = 0;
    else
        x_label = assign_noise(x, x_label');
        VRC(num) = compute_CH_index(x, x_label');
    end
    threshold(num) = ei;
    num = num + 1;
    %break;
end

%     figure
%     hold all;
%     stem(VRC)

[V, I] = max(VRC);
ei = threshold(I);

clusterId = 0;
x_label = zeros(data_size, 1);
for i = order
    if RD(1, i) > ei
        if CD(1, i) <= ei
            clusterId = clusterId + 1;
            x_label(i) = clusterId;
        else
            x_label(i) = 0;
        end
    else
        x_label(i) = clusterId;
    end
end

x_label = assign_noise(x, x_label');
x_label = x_label';

%     figure
%     hold all;
%     plot([1 : 1 : X_range], ones(1, X_range).*ei, 'r', 'LineWidth', 1.5);
%     %stem(CD(1, order), 'r', 'LineWidth', 1.5);
%     stem(RD(1, order), 'b', 'LineWidth', 1.5);
%     grid;
%     box on;
%     set(gca,'FontName','Times New Roman');
%     set(gca,'FontSize', 20);
%     set(gca, 'FontWeight', 'bold');
%     set(gca,'XLim',[1 X_range]);
%     set(gca,'YLim',[0 Y_range]);
%     xlabel('Object Index');
%     ylabel('Reachability Distance');
end





