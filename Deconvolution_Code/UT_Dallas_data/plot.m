clc; close all; clear all;

load('s_cvxeda2.mat');

sub_fron_neuro = [1, 5, 8, 9, 12, 16]; %Frontiers in Neuroscience subjects

for i = sub_fron_neuro
    figure,
    plot(s(i).x);
end