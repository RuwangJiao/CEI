function CEI = Cal_CEI(offSubpop, opt, Succ_best_sumv, num_equal, num_inequal, centers, best_feasi_flag)
    global DB
    global integral_data
    global Dmodel
    global index_dist
    % 根据子代个体距离各聚类中心的距离排序来判断使用哪个局部近似模型
    Tem_distance    = sum((repmat(offSubpop, size(centers, 1), 1) - centers).^2, 2);
    Distance        = Tem_distance';
    [~, index_dist] = min(Distance); % index_dist(1,:) 就是对应的第几个聚类中心的索引值，表明针对 offSubpop(i,:)个体，应该用第 index_dist(1,:) 个局部近似模型来预筛
    switch best_feasi_flag
        case 0
            g_min             = min(DB.G_value);
            integral_data.x   = offSubpop;
            integral_data.Num = 2*num_equal + num_inequal;
            term1             = integral(@(z) Fai(z), 0, g_min);
            improv_y          = 0;
            PI_total          = 1;
            for j = 2:(2*num_equal + num_inequal + 1)
                probImp1(j-1) = probImp_dace(offSubpop, improv_y, Dmodel(index_dist, j).dmodel);
                PI_total      = PI_total*probImp1(j-1);
            end
            term2 = g_min*PI_total;
            CEI   = term1 - term2;
        case 1
            improv_y = 0;
            ExpImp1  = ExpImp_dace(offSubpop, opt.y, Dmodel(index_dist, 1).dmodel);
            PI_total = 1;
            for j    = 2:(2*num_equal + num_inequal + 1)
                probImp1(j-1) = probImp_dace(offSubpop, improv_y, Dmodel(index_dist, j).dmodel);
                PI_total = PI_total*probImp1(j-1);
            end
            CEI = ExpImp1*PI_total;
    end
end