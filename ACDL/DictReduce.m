function [ D_new, T_new ] = DictReduce( D, T, param )
%%   direct reduction
%   Input:
%       D - dictionary, m-by-n array, n atoms and m dimension 
%       T - label
%            
%   Outputs:  
%       D_new - new dictionary
%       T_new - new label
%       param.dist - parameter for distance reduction
%       param.class - parameter for distinction reduction
%
%   Author: Yang Song (ysong18@utk.edu)
%
%% check input
if nargin < 2
    error('Not enough input arguments!')
elseif nargin == 2
    param.dist = .4;
    param.class = .5;
else
    if ~isfield(param, 'dist')
        param.dist = .1;
    end
    if ~isfield(param, 'class')
        param.class = .5;
    end
end
if size(D, 2) ~= size(T, 2)
 	error('Dismatching: D and T')
end

%% compute convex hull
% [K, V] = convhulln(D');
[m, n] = size(D);
% ind = unique(K(:));
% V = max(max(pdist2(D', D')));
ind = 1:n;
D_new = [];
T_new = [];

%% check correlation and distinction within class
D_temp = D(:, ind);
T_temp = T(:, ind);
[~, labels] = max(T_temp);
for i=1:length(unique(labels)) % vertices in a class
    tempT = T_temp(:, labels==i);
    tempD = D_temp(:, labels==i);
    if size(tempD, 2) < 2  % single vertex
        D_new = [D_new , tempD];
        T_new = [T_new , tempT];
        continue    
    end
    ind_kept = true(size(tempD,2), 1); % indicator of kept atoms
    % compute pair-wise distance
    % dist = pdist2(tempD', tempD') / V;
    dist = pdist2(tempD', tempD');
    dist_max = max( max(dist) );
    dist = dist / dist_max;
    [row, col] = find(dist>0 & dist<param.dist);
    if ~isempty(row)
        for j=1:length(row)
            gap = sort(tempT(:, [row(j), col(j)]));
            gap = gap(end, :) - gap(end-1, :);
            if gap(1) < gap(2)
                ind_kept(row(j)) = 0;
            else
                ind_kept(col(j)) = 0;
            end
        end
    end
    % compute distinctions 
    temp = sort(tempT);
    d = temp(end,:) - temp(end-1,:);
    ratio = d ./ max(d);
    ind_kept(ratio<param.class) = 0;
    D_new = [D_new , tempD(:, ind_kept)];
    T_new = [T_new , tempT(:, ind_kept)];
end

