function DB = update_database(xnew, addg, addh, f, addfit, num_inequal, delta, DB)
    %% Update database
    DB.Total_x       = [DB.Total_x; xnew];
    DB.Total_fit     = [DB.Total_fit; addfit]; 
    DB.Total_inequal = [DB.Total_inequal; addg];
    DB.Total_equal   = [DB.Total_equal; addh]; 
    DB.Total_obj     = [DB.Total_obj; f];
    
    % remove duplicated points
    Dis_eps       = 5e-8;
    index_delete  = [];
    xx_ind        = DB.Total_x;
    nn            = size(xx_ind, 1);
    traindisttem  = pdist(xx_ind, 'euclidean');
    traindist_ind = squareform(traindisttem);
    traindist_ind = traindist_ind + eye(nn);
    m111          = traindist_ind <= Dis_eps;
    for ij = 2:size(m111, 1)
        delete_flag = 0;
        for jk = 1:ij - 1
            if m111(ij, jk) >= 1 
                index_delete = [index_delete; ij];
                delete_flag = 1;
                if delete_flag == 1
                    break;
                end
            end
        end
    end
    DB.Total_x(index_delete, :)   = [];
    DB.Total_fit(index_delete, :) = [];
    if isempty(DB.Total_inequal)  ~= 1
        DB.Total_inequal(index_delete, :) = [];
    end
    if isempty(DB.Total_equal) ~= 1
        DB.Total_equal(index_delete, :) = [];
    end
    DB.Total_obj (index_delete, :) = [];

    DB.vio = sum(max(DB.Total_fit(:, 2:end), 0), 2);
   
    if num_inequal< (size(DB.fit, 2) - 1)
        DB.G_add = max(DB.Total_fit(:, 2:end), 0);
    else
        DB.G_add = max(DB.Total_fit(:, 2:end), 0);
    end        
    DB.G_value = max(DB.G_add,[],2); 
end