%% A Demo of Automatic Compact Dictionary Learning (ACDL)
% written by Yang Song, ysong18@utk.edu
% [CITATION]:
% Y. Song, Z. Zhang, L. Liu, A. Rahimpour, and H. Qi
% Dictionary Reduction: Automatic Compact Dictionary Learning for Classification
% Asian Conference on Computer Vision (ACCV), 2016.
% 
clc; clear; close all

%% add path of algorithm and data
addpath('ACDL');
addpath('data')

%% load data
load random_samples; % Y: 2-D observations, G: one-hot labels (two classes)
rng(0, 'v5uniform'); % reset random sampler

%% ACDL
param.Y = Y; 
param.G = G;
param.Dinit = Y; % initial dictionary 
[D, X, Yn, Err, T] = ACDL(param);

%% plot
param.D = D;
param.T = T;
param.title = 'ACDL';
plotDict(param);
