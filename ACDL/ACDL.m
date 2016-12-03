function [ D, X, Yn, Err, T ] = ACDL( param )
%   Automatic Compact Dictionary Learning (ACDL) for classification
%   objective function: 
%           min ||Y-DX||_f + gamma||G-TX ||_f + lambda*||X||_1 
%           s.t. sum||X||=1, X>=0
%   Input:
%       param.Y - observation matrix (M-by-N). Each column is an observation  
%       param.lambda - sparse penalty parameter
%       param.itermax - maximum iternation for searching optimal solution
%       param.threshold - Forbenis norm error (reconstruction error)
%       param.N_dict - number of atoms in the dictionary
%       param.Dinit - initial dictionary
%       param.Xinit - initial sparse coefficient
%       param.sum2one - flag of sum-to-one constrain. 1-enable (default), 0-disable
%       param.nonneg - flag of non-negative constrain. 1-enable (default), 0-disable
%       param.normcol - flag of dictionary normalization. 1-enable, 0-disable (default)
%       param.phi -  effect of sum-to-one constrain, default=1
%       param.gamma - classification error penalty parameter, default=1
%       param.G - label of observations (N columns). 
%       param.T - label of dictionary 
%       param.dist - from 0 to 1, normalized distance between two points. A
%       distance under this value indicated strong similarity.
%            
%   Outputs:  
%       D - leared dictionary
%       X - leared sparse coeffient
%       Yn - Yn=D*X, the reconstructed obserations
%       Err - reconstruction error, the forbenis norm between Yn and Y
%       T - estimated label of dictionary
%       
%
%   Author: Yang Song (ysong18@utk.edu)
%   Citation:
%   Y. Song, Z. Zhang, L. Liu, A. Rahimpour, and H. Qi
%   Dictionary Reduction: Automatic Compact Dictionary Learning for Classification
%   Asian Conference on Computer Vision (ACCV), 2016.
%  
%% check input arguments
if nargin < 1
    error('Not enough input arguments!'); end
if ~isfield(param, 'gamma')
    param.gamma = 1; end
if ~isfield(param, 'G')
    error('Missing param.G'); end
if ~isfield(param, 'Y')
    error('Missing param.Y'); 
else
    if size(param.Y, 2) ~= size(param.G, 2)
        error('Dimension dismatching: param.Y and param.G'); end
    paramDL.Y = [ param.Y ; param.gamma*param.G ]; 
end
if isfield(param, 'Dinit')
    paramDL.N_dict = size(param.Dinit, 2);
    paramDL.Dinit = [param.Dinit ; param.gamma*param.G];
elseif isfield(param, 'N_dict')
    paramDL.N_dict = param.N_dict;
else
    error('Missing param.N_dict');   
end
if isfield(param, 'lambda')
    paramDL.lambda = param.lambda;
else
    paramDL.lambda = .1; 
end
if isfield(param, 'itermax')
    paramDL.itermax = param.itermax; end
if isfield(param, 'threshold')
    paramDL.threshold = param.threshold; end
if isfield(param, 'Xinit')
    paramDL.Xinit = param.Xinit; end
if isfield(param, 'sum2one')
    paramDL.sum2one = param.sum2one; end
if isfield(param, 'nonneg')
    paramDL.nonneg = param.nonneg; end
if isfield(param,'normcol')
    paramDL.normcol = param.normcol; end
if isfield(param, 'phi')
    paramDL.phi = param.phi; end
if ~isfield(param, 'dist')
    param.dist = .4; 
end

%% automatic compact learn dictionary
i_D_end = size(param.Y, 1);
i_T_start = i_D_end+1;
while 1
    % dictionary learning
    [~, labels] = max(paramDL.Dinit(i_T_start:end, :), [], 1);
    n_labels = length(unique(labels));
    [D, X] = DirectDL(paramDL);
    D_new = D(1:i_D_end, :);
    T_new = D(i_T_start:end, :);
    [~, labels_new] = max(T_new, [], 1);
    n_labels_new = length(unique(labels_new));
    if n_labels ~= n_labels_new
        for i=1:n_labels
            if sum(labels_new==i) == 0
                indices = find(labels==i);
                D_new = [ D_new , paramDL.Dinit(1:i_D_end, indices) ];
                T_new = [ T_new , paramDL.Dinit(i_T_start:end, indices) ];
            end
        end
    end
    % dictionary reduction
    fprintf('\nDictsize of D input:%i\n',size(D_new,2));
    [D_out, T_out] = DictReduce(D_new, T_new, param);
    fprintf('\nDictsize of D output:%i\n',size(D_out,2));
    change = abs(size(D_out, 2)-size(D_new, 2));
    if (change < 3) || (size(D_new,2) < size(param.G,1));
        break; 
    else
        % option 1: update dictionary
        paramDL.Dinit = [D_out ; T_out];
        % option 2: update number of endmembers
        % paramDL.N_dict = size(D_out, 2);
    end
end

D = D_new;
T = T_new ./ param.gamma;
Yn = D * X;
Err = norm( (param.Y - Yn), 'fro');     

