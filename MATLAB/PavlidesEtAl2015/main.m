function main

% for more information see readme.txt

% load data
var = load('weights.mat');

[~, Features_opt] = minfunction(var.weights(1,:));
generate_fig(Features_opt, 1)

[~, Features_opt] = minfunction(var.weights(2,:));
generate_fig(Features_opt, 2)


