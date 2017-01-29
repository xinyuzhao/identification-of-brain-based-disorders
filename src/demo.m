clear all;
close all;
clc;

input_f = load('../data/connectivity_features/PTSD_connectivity.mat');
input_data = input_f.connectivities;

reference_label = load('sample_label.csv');

[bsetScore, bestLabel, bestFeature] = runPipeline(1, input_data, reference_label);
