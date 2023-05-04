function DB = build_database(initial_p, g, h, f, fit,num_inequal, delta)
    DB.x               = initial_p;
    DB.inequal         = g;
    DB.equal           = h;  % 可能不存在等式约束  如果 h=[] 表示不存在等式约束
    DB.obj             = f;
    DB.fit             = fit;
    DB.Total_x         = initial_p;
    DB.Total_inequal   = g;
    DB.Total_equal     = h;  % 可能不存在等式约束  如果 h=[] 表示不存在等式约束
    DB.Total_obj       = f;
    DB.Total_fit       = fit;
    DB.succ_delete_num = 0;
    DB.num_equal       = size(h, 2);
    DB.num_inequal     = size(g, 2);
    DB.vio             = sum(max(DB.fit(:, 2:2 + num_inequal-1), 0), 2);
    if num_inequal < (size(DB.fit,2) - 1)
        DB.vio = DB.vio + sum(max(abs(DB.fit(:, 2 + num_inequal:end)) - delta, 0), 2);
    end
    
    if num_inequal < (size(DB.fit, 2) - 1)
        DB.G_add = [max(DB.fit(:, 2:2 + num_inequal - 1), 0)  max(abs(DB.fit(:, num_inequal + 2:end)) - delta, 0)];
    else
        DB.G_add = max(DB.fit(:, 2:end), 0);
    end
            
    DB.G_value = max(DB.G_add, [], 2);       

    if isempty(h) == 1
        DB.exist_equal = 1; 
    else
        DB.exist_equal = 0;
    end
    
    if isempty(g) == 1
        DB.exist_inequal = 1;
    else
        DB.exist_inequal = 0;
    end
end