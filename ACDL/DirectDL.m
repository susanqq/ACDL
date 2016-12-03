function[D, X, Yn, Err] = DirectDL( param )
%%   direct dictionary learning
%   objective function: min ||Y-DX||_f + lambda*||X||_1 
%                       s.t. sum||X||=1, X>=0
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
%       param.phi -  effect of sum-to-one constrain. default=1.
%            
%   Outputs:  
%       D - leared dictionary
%       X - leared sparse coeffient
%       Yn - Yn=D*X, the reconstructed obserations
%       Err - reconstruction error, the forbenis norm between Yn and Y
%       
%
%   Author: Yang Song (ysong18@utk.edu)
%
%% check input arguments
if nargin < 1
    error('Not enough input arguments!'); end
if ~isfield(param, 'Y')
    error('Missing param.Y'); end
n_rows = size(param.Y, 1);
if isfield(param, 'Dinit') 
    param.N_dict = size(param.Dinit, 2);
elseif ~isfield(param, 'N_dict')
    error('Missing param.N_dict');
else
    D0 = rand(n_rows, param.N_dict); 
    % normalization:  divide each column by column norm
    param.Dinit = D0 ./ (ones(n_rows,1)*sqrt(sum(D0.^2)));
    param.N_dict = size(param.Dinit, 2);
end
if ~isfield(param, 'lambda') 
    error('Missing param.lambda'); end
if ~isfield(param, 'itermax')
    param.itermax=50; end
if ~isfield(param, 'threshold')
    param.threshold=1e-3; end
if ~isfield(param, 'Xinit')
    param.Xinit = rand(param.N_dict, size(param.Y,2) ); end 
if ~isfield(param, 'sum2one')
    param.sum2one = 1; end
if ~isfield(param, 'nonneg')
    param.nonneg =1; end
if ~isfield(param,'normcol')
    param.normcol =0; end
if ~isfield(param, 'phi')
    param.phi = 1; end

if size(param.Dinit, 2) ~= size(param.Xinit, 1)
    error('Dimension dismatching: param.Dinit and param.Xinit.'); end
if size(param.Dinit, 1) ~= size(param.Y, 1) || size(param.Xinit, 2) ~= size(param.Y, 2)
    error('Dimension dismatching: Y and param.Dinit or .Xinit'); end

%% add sum-to-one constrain to Y and D
Y = param.Y;
D = param.Dinit;
X = param.Xinit;
if param.sum2one 
    tempY = ones(1, size(Y,2)) * param.phi;
    tempD = ones(1, size(D,2)) * param.phi;
    Y = [Y ; tempY];
    D = [D ; tempD];  
end
    
%% update dictionary D and sparse coefficient X
for i=1:param.itermax
    % sum-to-one constrain
    if param.sum2one
        D(end,:) = tempD;
        Y(end,:) = tempY;
    end
    % stepsize parameter for gradient descent
    tau1 = 1/norm(D*D');
    tau2 = 1/norm(X*X');
    Xo = X;
    Do = D;
    DoXo = Do * Xo;
    gradX=Do' * ( Y - DoXo );
    gradD = ( Y - DoXo ) * Xo';
    X = Xo + tau1 * gradX;
    X= prox_l1(X, param.lambda * tau1);
    if param.nonneg 
        X = max(X, 0); end
    D = Do + tau2 * gradD;
    % make all atom in unit-norm D
    if param.normcol && param.sum2one
        D = D(1:end-1, :);
        Y = Y(1:end-1, :);
        temp = max(1, sqrt(sum(D .^ 2)));
        index = find(temp>1);
        D(:, index) = D(:,index)./(ones(size(D,1),1)*temp(index));
        D=[D;tempD];  
        Y=[Y;tempY];
    elseif param.normcol
        temp=max(1,sqrt(sum(D.^2)));
        index=find(temp>1);
        D(:,index) = D(:,index)./(ones(size(D,1),1)*temp(index));
    end
     % stop rule
     error_D=norm(Do - D,'fro') / norm(Do,'fro') ;
    if error_D< param.threshold
        break
    end
    fprintf('\niteration:%f,error:%f',i,error_D);
end % end of iteration
% remove sum-to-one constrain
if param.sum2one
    Y(end, :) = [];
    D(end, :) = [];
end
Yn = D * X;
[m, n] = size(Y);
Err=norm( (Y-Yn), 'fro') / (m*n);     

