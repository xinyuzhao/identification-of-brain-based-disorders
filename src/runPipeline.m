% -------------------------------------------------------------------------
% Function:
% [bestScore, bestLabel, bestFeature]=runPipeline(method_name, input_data, reference_label)
% -------------------------------------------------------------------------
% Target:
% Determine the relevent subset of features and cluster samples into
% different groups simultaneously.
% -------------------------------------------------------------------------
% Input:
% method_name - selected clustering method; 0: hierarchial clustering,
% 1: DPC, 2: OPTICS
% input_data - m by n data matrix; m: number of samples, n: number of
% features
% reference_label - m by 1 vector; Optional input paramter, default is []
% -------------------------------------------------------------------------
% Output:
% bestScore - float number; best score
% bestLabel - m by 1 vector
% bestFeature - vector specifying the selected features indices
% -------------------------------------------------------------------------
% References:
% [1] Dy, J.G., Brodley, C.E., 2004. Feature Selection for Unsupervised
% Learning. J. Mach. Learn. Res. 5., 845-889.
% -------------------------------------------------------------------------
% Written by Xinyu Zhao
% AU MRI research center, Auburn University
% Janunary, 2017
% https://github.com/xinyuzhao/identification-of-brain-based-disorders

function [bestScore, bestLabel, bestFeature] = ...
    runPipeline(method_name, input_data, reference_label)
if nargin < 3
    reference_label = [];
end

numFeature = size(input_data, 2);
%% Parameter Setting
% MaxNumClust: maximum number of clusters
% MaxIteration: maximum number of iterations
% MaxPercentFeature: maximum percentage of features that can be selected
% popSize: population size
% pm: permutaion rate
% goodSolution: percentage of top solutions that will be selected into next
% generation
% givenPotentialSolution: if true load from file, else randomly generate
MaxNumClust = 5;
MaxIteration = 1000;
MaxPercentFeature = 0.5;
popSize = 400;
pm = .3;
goodSolution = .2;
givenPotentialSolution = false;
rng('shuffle')

%% Initialize
pop = InitialPop(numFeature, popSize, MaxPercentFeature, givenPotentialSolution);

for iter = 1 : MaxIteration
    %% Crossover
    kids = Crossover(pop, popSize, numFeature);
    
    %% Mutation
    kids = Mutation(kids, popSize / 2, numFeature, pm);
    
    %% Evaluation
    totalPop = [pop; kids];
    save('result/totalPop.mat', 'totalPop');
    [bestScore, bestLabel, bestFeature, sorted_sims] = Evaluation(input_data, reference_label, MaxNumClust, totalPop, iter, method_name);
    
    %% SelectPop
    pop = SelectPop(totalPop, sorted_sims, popSize, numFeature, goodSolution);
end
end

function pop = InitialPop(chromosomeLength, popSize, MaxPercentFeature, givenPotentialSolution)
if givenPotentialSolution == true
    input = load('result/latest_pop.mat');
    pop = input.newPop;
else
    pop = zeros(popSize, chromosomeLength);
    for i = 1 : popSize
        while(nnz(pop(i, :))==0)
            for j = 1 : chromosomeLength
                if (rand(1) < MaxPercentFeature)
                    pop(i, j) = 1;
                end
            end
        end
    end
end
end

function kids = Crossover(pop, popSize, numFeature)
parentIndex = randperm(popSize);
kids = zeros(popSize, numFeature);
for i = 1 : 2 : popSize - 1
    p1 = pop(parentIndex(i), :);
    p2 = pop(parentIndex(i + 1), :);
    V1 = randi([1 numFeature]);
    kids(i, :) = [p1(1 : V1), p2(V1 + 1 : end)];
    kids(i + 1, :) = [p2(1 : V1), p1(V1 + 1 : end)];
end
end

function kids = Mutation(kids, kidsSize, numFeature, pm)
p1 = pm;
p0 = pm;
for i = 1 : kidsSize
    for j = 1 : numFeature
        r = rand(1);
        if kids(i, j) == 1 && r < p1 && nnz(kids(i, :)) ~= 1
            kids(i, j) = 0;
        elseif kids(i, j) == 0 && r < p0
            kids(i, j) = 1;
        end
    end
end
end

function [bestScore, bestLabel, bestFeature, sorted_sims] = Evaluation(data, reference, MaxNumClust, totalPop, iter, method_name)
totalNum = size(totalPop, 1);
sims = zeros(totalNum, 1);
bestScore = -1;

for i = 1 : totalNum
    picked_feature = find(totalPop(i, :) == 1);
    if nnz(totalPop(i, :)) ~= 0
        if method_name == 0
            [label, score] = cluster_hierarchical(data(:, picked_feature), MaxNumClust);
        elseif method_name == 1
            [label, score] = cluster_dp(data(:, picked_feature), MaxNumClust);
        elseif method_name == 2
            [label, score] = cluster_optics(data(:, picked_feature), MaxNumClust);
        else
            disp('Invalid method name');
            return
        end
        if ~isempty(reference)
            sims(i, 1) = similarity(reference, label);
        else
            sims(i, 1) = score;
        end
    else
        sims(i, 1) = 0;
    end
    if sims(i, 1) > bestScore
        bestScore = sims(i, 1);
        bestLabel = label;
        bestFeature = picked_feature;
    end
end

combine_sims = [[1 : totalNum]', sims];
sorted_sims = sortrows(combine_sims, -2);
numF = nnz(bestFeature);
disp([iter, bestScore, numF]);
save('result/latest_best.mat', 'bestFeature', 'bestScore', 'bestLabel');
end

function newPop = SelectPop(totalPop, sorted_sims, popSize, numFeature, goodSolution)
newPop = zeros(popSize, numFeature);
numSolution = size(totalPop, 1);
selectedSolution = sorted_sims([1 : popSize * goodSolution, ...
    numSolution - popSize * (1 - goodSolution) + 1 : numSolution], 1);
newPop(:, :) = totalPop(selectedSolution, :);
save('result/latest_pop.mat', 'newPop');
end

