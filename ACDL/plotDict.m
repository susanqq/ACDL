function [ f ] = plotDict( param )
%%  plot dictionary, observation, reconstruction, etc. in a 2D space
%
% param. Y - observation
% param. D - dictionary
% param. Yn - reconstruction
% param. T - labels of dictionary
% param. G - labels of observation 
% param.title - title of the figure
%
% Auther: Yang Song (ysong18@utk.edu)
% Created on 8/6/2015
%

%%  check input arguments
if nargin < 1
    error('Not enough input arguments!'); end
if ~isfield(param, 'Y')
    error('Missing param.Y'); end
if ~isfield(param, 'D')
    error('Missing param.D'); end
if isfield(param, 'Yn')
    flag_Yn = 1;
    if size(param.Y, 2) ~= size(param.Yn, 2)
        error('Dismatching: param.Y and .Yn'); end
else
    flag_Yn = 0;
end
if isfield(param, 'T')
    flag_Dcluster = 1;
    if size(param.D, 2) ~= size(param.T, 2)
        error('Dismatching: param.D and .T'); end
else
    flag_Dcluster = 0;
end
if isfield(param, 'G')
    flag_Ycluster = 1;
    if size(param.Y, 2) ~= size(param.G, 2)
        error('Dismatching: param.Y and .G'); end
else
    flag_Ycluster = 0;
end


%% get painting materials
f = figure; hold on; 
colors = get(gca, 'colororder');
colors_labels = ['c', 'm']';
markers = {'o', 's', '^', 'x', '+', '*', 'd'};
lineWidth = 2;
circleArea = 50;


%% plot reconstruction error
Y = param.Y;
if flag_Yn
    Yn = param.Yn;
    for i=1:size(Y, 2)
        plot( [Y(1,i), Yn(1,i)],  [Y(2,i), Yn(2,i)], 'k:' )
    end
end

%% plot observation Y
if flag_Ycluster
    G = param.G;
    for i=1:size(G, 2)
        [~, label] = max( G(:, i) ); 
        scatter( Y(1, i), Y(2, i), circleArea, markers{label},'MarkerEdgeColor', 'k', ...
              'MarkerFaceColor', colors(label, :) )
    end
else
    scatter( Y(1, :), Y(2, :), circleArea, markers{1}, 'MarkerEdgeColor', 'k', ...
              'MarkerFaceColor', colors(1, :) )
end

%% plot reconstruction Yn
if flag_Yn
    Yn = param.Yn;
    if flag_Ycluster
        G = param.G;
        for i=1:size(G, 2)
            [~, label] = max( G(:, i) ); 
            scatter( Yn(1, i), Yn(2, i), circleArea, markers{label}, ...
                'MarkerEdgeColor', colors(label, :), ...
                'linewidth', lineWidth)
        end
    else
        scatter( Yn(1, :), Yn(2, :), circleArea, markers{1}, ...
            'MarkerEdgeColor', colors(1, :), ...
            'linewidth', lineWidth)
    end
end

%% plot convex hull
D = param.D(1:2, :);
dn = size(D, 2);
if dn > 2
    k = convhull(D(1,:)', D(2,:)');
else
    k = 1:dn;
end
plot( D(1,k), D(2,k), 'k--', 'linewidth', lineWidth );

%% plot dictionary D
if flag_Dcluster
    T = param.T;
    for i=1:size(T, 2)
        [~, label] = max( T(:, i) ); 
        scatter( D(1, i), D(2, i), circleArea*10, 'h', 'MarkerEdgeColor', 'k', ...
              'MarkerFaceColor', colors_labels(label, :) )
    end
else
    scatter( D(1, :), D(2, :), circleArea*10, 'h', 'MarkerEdgeColor', 'k', ...
              'MarkerFaceColor', colors(2, :) )
end

%% font size
fontSize_L = 20;
fontSize_S = 18;
set(gca, 'fontsize', fontSize_S);
if isfield(param, 'title')
    title(param.title, 'fontsize', fontSize_L);
end
grid on

end % end of function










