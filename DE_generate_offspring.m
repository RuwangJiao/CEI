function [population, bestvalue, bestind] = DE_generate_offspring(NP, lu, n, F, CR, opt, Succ_best_sumv, num_equal, num_inequal, centers, best_feasi_flag, Maxt)
    %% DE Parameters
    lb    = lu(1, :);
    ub    = lu(2, :);
    MaxIt = Maxt;      % Maximum Number of Iterations最大迭代次数
    nPop  = NP;        % Population Size种群大小

    %% Initialization
    empty_individual.Position = [];
    empty_individual.Cost     = [];
    BestSol.Cost              = inf;
    pop                       = repmat(empty_individual, nPop, 1);
    for i = 1:nPop
         pop(i).Position = rand(1, size(lu,2)).*(lu(2,:) - lu(1,:)) + lu(1,:);
         pop(i).Cost     = Cal_CEI(pop(i).Position, opt, Succ_best_sumv, num_equal, num_inequal, centers, best_feasi_flag);
         pop(i).Cost     = -pop(i).Cost;
        if pop(i).Cost < BestSol.Cost
            BestSol.Cost = pop(i).Cost; 
            BestSol.ind  = pop(i).Position;
        end
    end
    BestCost = zeros(MaxIt, 1);
    %% DE Main Loop
    beta = F;
    pCR  = CR;
    FES  = nPop;
    for it = 1:MaxIt
        for i = 1:nPop
            x       = pop(i).Position;     
            A       = randperm(nPop);
            A(A==i) = [];
            a       = A(1);
            b       = A(2);
            c       = A(3);
            y       = pop(a).Position + beta.*(pop(b).Position - pop(c).Position);
            VarMin  = repmat(lb, size(y, 1), 1);
            VarMax  = repmat(ub, size(y, 1), 1);
            y       = max(y, VarMin);
            y       = min(y, VarMax);   
            % Crossover
            z       = zeros(size(x));
            j0      = randi([1 numel(x)]);
            for j = 1:numel(x)
                if j == j0 || rand <= pCR
                    z(j) = y(j);
                else
                    z(j) = x(j);
                end
            end

             NewSol.Position = z;
             NewSol.Cost     = Cal_CEI(NewSol.Position, opt, Succ_best_sumv, num_equal, num_inequal, centers, best_feasi_flag);
             NewSol.Cost     = -NewSol.Cost;
             FES             = FES + 1;

            if NewSol.Cost < pop(i).Cost
                pop(i) = NewSol;
                if pop(i).Cost < BestSol.Cost
                   BestSol.Cost = pop(i).Cost;
                   BestSol.ind  = pop(i).Position;
                end
            end
        end

        % Update Best Cost更新最佳适应度值
        BestCost(it) = BestSol.Cost;
        bestvalue    = BestSol.Cost;
        bestind      = BestSol.ind;

        for i = 1:nPop
            population(i, :) = pop(i).Position;
        end
    end
end