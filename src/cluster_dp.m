% -------------------------------------------------------------------------
% Function:
% [cl, V] = cluster_dp(data, max_num)
% -------------------------------------------------------------------------
% Target:
% Cluster data using density peak clustering method
% -------------------------------------------------------------------------
% Input:
% data - m by n data matrix; m: number of samples, n: number of
% max_num - integer; maximum number of clusters
% -------------------------------------------------------------------------
% Output:
% cl - m by 1 vector
% V - float number; clustering score
% -------------------------------------------------------------------------
% References: 
% [1] Rodrugyezm A., Laio, A., 2014. Clustering by fast search and find of
% density peaks. Science 344, 1492-6.
% -------------------------------------------------------------------------
% Written by Xinyu Zhao
% AU MRI research center, Auburn University
% Janunary, 2017
% https://github.com/xinyuzhao/identification-of-brain-based-disorders

function [cl, V] = cluster_dp(data, max_num)
num_subject = size(data, 1);
xx = pdist(data)';

ND = num_subject;
NL = num_subject;

N = size(xx, 1);

percent=2.0;

position=round(N * percent/100);
sda=sort(xx);
dc=sda(position);
if dc == 0
    loc = min(find(sda ~= 0));
    dc = sda(loc);
end

rho = zeros(ND, 1);

% Gaussian kernel

for i = 1 : ND - 1
    for j = i + 1 : ND
        distance = xx((i - 1)*(num_subject - i/2) + j - i);
        rho(i)=rho(i)+exp(-(distance/dc)*(distance/dc));
        rho(j)=rho(j)+exp(-(distance/dc)*(distance/dc));
    end
end

maxd=max(xx);

[rho_sorted, ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;

for ii = 2 : ND
    delta(ordrho(ii)) = maxd;
    for jj = 1 : ii - 1
        tmp_i = ordrho(ii);
        tmp_j = ordrho(jj);
        if tmp_i > tmp_j
            tmp = tmp_i;
            tmp_i = tmp_j;
            tmp_j = tmp;
        end
        if tmp_i == tmp_j
            distance = 0;
        else
            distance = xx((tmp_i - 1)*(num_subject - tmp_i/2) + tmp_j - tmp_i);
        end
        if(distance < delta(ordrho(ii)))
            delta(ordrho(ii)) = distance;
            nneigh(ordrho(ii)) = ordrho(jj);
        end
    end
end
delta(ordrho(1))=max(delta(:));

loc = find(nneigh == 0);
    
% mean(rho)    
rhomin_list = mean(rho) : 0.5 : min(mean(rho) + 10, max(rho) * 0.8);
    deltamin_list = max(delta) * 0.99 : -mean(delta) * 0.5 : mean(delta) * 0.5;
    
    opt_solution = zeros(1, 3);
    index = 1;
    for num_c = 2 : max_num
        for rhomin = rhomin_list
            for deltamin = deltamin_list
                NCLUST=0;
                cl = ones(1, ND) * -1;
                
                for i = 1 : ND
                    if ( (rho(i) > rhomin) && (delta(i) > deltamin))
                        NCLUST = NCLUST+1;
                        cl(i) = NCLUST;
                        icl(NCLUST)=i;
                    end
                end
                
                %assignation
                for i = 1 : ND
                    if (cl(ordrho(i))==-1)
                        cl(ordrho(i))=cl(nneigh(ordrho(i)));
                    end
                end
                
                if (NCLUST > 1)
                    bord_rho = zeros(1, ND);
                    for i=1:ND-1
                        for j=i+1:ND
                            distance = xx((i - 1)*(num_subject - i/2) + j - i);
                            if ((cl(i)~=cl(j))&& (distance <= dc))
                                rho_aver=(rho(i)+rho(j))/2.;
                                if (rho_aver>bord_rho(cl(i)))
                                    bord_rho(cl(i))=rho_aver;
                                end
                                if (rho_aver>bord_rho(cl(j)))
                                    bord_rho(cl(j))=rho_aver;
                                end
                            end
                        end
                    end
                end %end of if (NCLUST > 1)
                if NCLUST == num_c
                    tmp = compute_CH_index(data, cl);
                    opt_solution(index, 1) = rhomin;
                    opt_solution(index, 2) = deltamin;
                    opt_solution(index, 3) = tmp;
                    index = index + 1;
                end
            end
        end
    end
    
    [V, I] = max(opt_solution(:, 3));
    rhomin = opt_solution(I, 1);
    deltamin = opt_solution(I, 2);
    NCLUST=0;
    cl = ones(1, ND) * -1;
    
    for i = 1 : ND
        if ( (rho(i) > rhomin) && (delta(i) > deltamin))
            NCLUST = NCLUST+1;
            cl(i) = NCLUST;
            icl(NCLUST)=i;
        end
    end
    
    %assignation
    for i = 1 : ND
        if (cl(ordrho(i))==-1)
            cl(ordrho(i))=cl(nneigh(ordrho(i)));
        end
    end
    
    if (NCLUST > 1)
        bord_rho = zeros(1, ND);
        for i=1:ND-1
            for j=i+1:ND
                distance = xx((i - 1)*(num_subject - i/2) + j - i);
                if ((cl(i)~=cl(j))&& (distance <= dc))
                    rho_aver=(rho(i)+rho(j))/2.;
                    if (rho_aver>bord_rho(cl(i)))
                        bord_rho(cl(i))=rho_aver;
                    end
                    if (rho_aver>bord_rho(cl(j)))
                        bord_rho(cl(j))=rho_aver;
                    end
                end
            end
        end
    end %end of if (NCLUST > 1)
    
%     figure
%     hold all;
%     plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
%     plot(rho(:), ones(size(rho, 1), 1).*deltamin, 'r')
%     title ('Decision Graph','FontSize',15.0)
%     xlabel ('\rho')
%     ylabel ('\delta')
    
    cl = cl';

end

