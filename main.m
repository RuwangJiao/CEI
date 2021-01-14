% Constrained Expected Improvement (CEI) 
%
% Reference: Jiao R, Zeng S, Li C, Jiang Y, Jin Y. A complete expected improvement criterion for Gaussian process 
%            assisted highly constrained expensive optimization[J]. Information Sciences, 2019, 471: 80-96.
%
% ruwangjiao@gmailcom  
% Last update January 14, 2021

clear all;
clc;
addpath(genpath('dace'));
global Dmodel
global DB
global integral_data

problem = 3;
fprintf('G%s problem \n', num2str(problem));
[n, lu, AA]     = cec2006Suite(problem);
initial_popsize = 11*n - 1;
initial_p       = repmat((lu(2,:) - lu(1,:)), initial_popsize, 1).*lhsdesign(initial_popsize, n, 'criterion','maximin','iterations',100) + repmat(lu(1,:), initial_popsize, 1);
mu_equal        = 0;  
delta_equal     = 1e-4; 
[fit,h1,g,f]    = fitness(initial_p, problem, delta_equal, AA);   
[h, num_inequal, num_equal,fit] = correct_obj_con_equ(h1, g, fit);
[~,~,~,Tem_fit] = correct_obj_con(h1, g, fit);

%% Training set
DB             = build_database(initial_p, g, h, f, fit, num_inequal, delta_equal);
FES            = initial_popsize;
Succ_best_sumv = [];
[opt, Succ_best_sumv, best_feasi_flag] = BestOne_2(initial_p, fit,num_inequal, num_equal, delta_equal, FES,Succ_best_sumv);
gen    = 1;
MaxFES = 50*n;
while FES <= MaxFES  
   %% fuzzy c-means
    L1              = 80;
    L2              = 20;
    data_c_means    = DB.Total_x;
    N               = size(data_c_means,1);
    Nc              = 1 + max(ceil((N-L1)/L2), 0);
    cluster_options = [2.0 NaN 0.05 0];
    [centers, U]    = fcm(data_c_means, Nc, cluster_options); 
    for i = 1:Nc
        [~, sort_index(i).index1] = sort(U(i, :), 'descend');
        if Nc <= 1
            index2(i).index1 = sort_index(i).index1(1:N);
        else
            index2(i).index1 = sort_index(i).index1(1:L1);
        end
    end
    
   %% Initialize hyper parameters 
    theta = 10*rand(1, n);
    lob   = (1e-4)*ones(1, n);
    upb   = 10*ones(1, n);
    Total_y_sample = [];
    Total_x_sample = [];
    numm  = (2*num_equal + num_inequal + 1);
    
    % Training models
    for i = 1:Nc
        X_sample(i).sample = DB.Total_x(index2(i).index1, :);
        Y_sample(i).sample = DB.Total_fit(index2(i).index1, :);
        YY_sample          = Y_sample(i).sample;
        XX_sample          = X_sample(i).sample;
        for j = 1:numm
            y_sample = YY_sample(:, j);
            x_sample = XX_sample;
            [dmodel, perf] = dacefit(x_sample, y_sample, @regpoly1, @corrgauss,theta, lob, upb);
            Dmodel(i,j).dmodel = dmodel; 
        end
    end
  
   %% DE generation offspring
    NP   = 30;
    F    = 0.5;
    CR   = 0.9;
    Maxt = 500;
    [offSubpop,~,bestind] = DE_generate_offspring(NP, lu, n, F, CR, opt, Succ_best_sumv, num_equal, num_inequal, centers, best_feasi_flag, Maxt);
  
    Succ_ind(gen, :)                       = bestind;
    [addfit, addh1, addg, addf]            = fitness(bestind, problem, delta_equal, AA);
    [addh, num_inequal, num_equal, addfit] = correct_obj_con_equ(addh1,addg,addfit);
    [~, ~, ~, Tem_addfit]                  = correct_obj_con(addh1, addg, addfit);
    DB                                     = update_database(bestind, addg, addh, addf, addfit, num_inequal, delta_equal, DB);
    FES = FES+1;
    [opt, Succ_best_sumv, best_feasi_flag] = BestOne_2(DB.Total_x, DB.Total_fit, num_inequal, num_equal, delta_equal, FES, Succ_best_sumv);  
    gen = gen + 1;
    fprintf('FEs=%d, BestObj=%f, BestVio=%f, PredObj=%f, PredVio=%f \n', FES-1, opt.y, min(DB.G_value),addfit(:,1), max(0,max(addfit(:,2:end))));
end
