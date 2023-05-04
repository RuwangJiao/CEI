function [opt, Succ_best_sumv, best_feasi_flag] = BestOne_2(x, y, num_inequal, num_equal, delta, FES, Succ_best_sumv)
    obj   = y(:, 1);
    vio   = sum(max(y(:, 2:end), 0), 2);
    index = find(vio == 0);
    if ~isempty(index)
        [opt.y, index2]           = min(obj(index));
        opt.x                     = x(index(index2,:), :);
        Succ_best_sumv.val(FES,:) = opt.y;
        Succ_best_sumv.fit(FES,:) = y(index(index2,:), :);
        Succ_best_sumv.ind(FES,:) = opt.x;
        best_feasi_flag           = 1;
    else
        opt.y                     = NaN;
        opt.x                     = NaN;
        [min_vio,index3]          = min(vio);
        Succ_best_sumv.val(FES,:) = y(index3, :);
        Succ_best_sumv.fit(FES,:) = y(index3, :);
        Succ_best_sumv.ind(FES,:) = x(index3, :);
        best_feasi_flag           = 0;
    end
end