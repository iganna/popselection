function obj_func = pytree_congruence(dist_x, freqs_x, dist_y, freqs_y, show_flag)
% pytree_congruence - This function defines the optimal Procrustes 
%          transformation to minimize the distance between X and Y
%          populatios represented by distance matrixes between sequences
%          and their frequences
%
% Inputs:
%    dist_x  - a distance matrix between sequences contained in X population
%    freqs_x - a vector of frequences of sequences in X population
%    dist_y  - a distance matrix between sequences contained in Y population
%    freqs_y - a vector of frequences of sequences in Y population
%    show_flag - true or false; to visualise the result of analysis or not
%
% Outputs:
%    obj_func - the optimal value between two gaussian mixture models
%
% Other m-files required: GaussMix2Distr_total,
%                         GaussMix2Distr_regul,
%                         procrustes_calc,
%                         cm_purple, cm_orange - two colormaps from
%                                               cm_space_plots.mat
%                         
% Author: Anna A. Igolkina
% email address: igolkinaanna11@gmail.com
% Last revision: 01-Jan-2018

%% Default parameters
sigma_scale = 1; % covariance matrix of Gaussian distributions
n_dim = 3;       % dimension of a space where X and Y population will be presented
d_scale_new = 3;    % The median distance between sequences after 
                 % the normalisation of distances between sequences
clear sigma
for i = 1:length(freqQuery)
    sigma(:,:,i) = diag(ones(n_dim,1))  * sigma_scale;
end
                 
%% Translation of the a population into the n_dim space

D_up = @(D) triu(D,1)
D_scale = @(D_up)median(D_up(D_up > 0))
D_scaled = @(D, D_scale, D_scale_new) D/D_scale * D_scale_new
D_norm = @(D, D_scale_new) D_scaled(D, D_scale(D_up(D)), D_scale_new)
MDS = @(D, D_scale_new, N_dim) mdscale(D_norm(D, D_scale_new), N_dim)


%% Apply translation
[X, stress] = MDS(dist_x, d_scale_new, n_dim);
sigma_x = sigma; 
fx = freqs_x/sum(freqs_x);

[Y_init, stress] = MDS(dist_y, d_scale_new, n_dim);
sigma_y = sigma; 
fy = freqs_y/sum(freqs_y);


% X = X(fx ~= 0, :);
% sigma_x = sigma_x(:,:,fx ~= 0);
% fx = fx(fx ~= 0);
% Y = Y(fy ~= 0, :);
% sigma_y = sigma_y(:,:,fy ~= 0);
% fy = fy(fy ~= 0);


%% Procrustes transformation

obj_func_set = [];
params_set = {};

for Y_it = {Y_init, -Y_init} % for mirror reflection
    Y = Y_it{1}
    my_obj = @(param)GaussMix2Distr_regul(param',X,Y,sigma_x, sigma_y, fx, fy, 0);
    problem = createOptimProblem('fmincon','x0', [rand(n_dim*(n_dim-1)/2,1) ; rand(n_dim ,1)*3] ,...
        'lb',[0 0 0 -Inf -Inf -Inf],...
        'ub', [2*pi 2*pi 2*pi Inf Inf Inf],...
        'objective',my_obj);

    ms = MultiStart;
    [params,~,~,~,~] = run(ms,problem,100);

    obj_func = GaussMix2Distr_total(params',X,Y,sigma_x, sigma_y, fx, fy)

    obj_func_set = [obj_func_set; obj_func];
    params_set = [params_set; {params}];

end

if obj_func_set(1) > obj_func_set(2)
    obj_func = obj_func_set(2);
    params = params_set{2};
else
	obj_func = obj_func_set(1);
    params = params_set{1};
end


%% Visualisation

if show_flag == true

    Z = procrustes_calc(Y, params'); 
    
    figure;
    hold on;
    colx = getcolormap(fx, cm_purple);
    for i = 1:size(X,1)
        plot3(X(i,1), X(i,2), X(i,3), 'o', 'Color', colx(i,:), 'MarkerSize', fx(i) * 100,  'LineWidth',2);
    end
    hold on; grid on;
    colz = getcolormap(fy, cm_orange);
    for i = 1:size(Z,1)
        plot3(Z(i,1), Z(i,2), Z(i,3), 'o', 'Color', colz(i,:), 'MarkerSize', fy(i) * 100, 'LineWidth',2);
    end

end
