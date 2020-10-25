clc; clear all; close all

% load bpopt libraries
addpath('bpopt')

% Solve example problem 4
prob = bp_create(5);  % create problem
prob.mu = 10;          % change the initial barrier term
sol = bpopt(prob);     % solve problem

