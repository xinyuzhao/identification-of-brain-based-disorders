% -------------------------------------------------------------------------
% Function:
% CH=compute_CH_index(data, label)
% -------------------------------------------------------------------------
% Target:
% Compute the Calinski-Harabasz index 
% -------------------------------------------------------------------------
% Input:
% data - m by n data matrix; m: number of samples, n: number of features
% label - m by 1 vector; 
% -------------------------------------------------------------------------
% Output:
% CH - float number;
% -------------------------------------------------------------------------
% References: 
% [1] Calinski, T., Harabasz, J., 1974. A Dendrite Method for cluster
% Analysis. Commun. Stat.-Simul. Comput. 3, 1-27.
% -------------------------------------------------------------------------
% Written by Guy Brys and Wai Yan Kong 
% May 2006
% http://wis.kuleuven.be/stat/robust.html

function CH = compute_CH_index(data, label)
[nrow, ncol] = size(data);

k = length(unique(label));
[st, sw, sb, cintra, cinter] = valid_sumsqures(data, label, k);

ssw = trace(sw);
ssb = trace(sb);

if k > 1
  CH = ssb/(k-1); 
else
  CH =ssb; 
end

CH = (nrow-k)*CH/ssw;
end