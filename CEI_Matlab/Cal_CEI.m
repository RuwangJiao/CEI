function CEI = Cal_CEI(offSubpop, opt, Succ_best_sumv, num_equal, num_inequal, centers, best_feasi_flag)
    global DB
    global integral_data
    global Dmodel
    global index_dist
    % �����Ӵ����������������ĵľ����������ж�ʹ���ĸ��ֲ�����ģ��
    Tem_distance    = sum((repmat(offSubpop, size(centers, 1), 1) - centers).^2, 2);
    Distance        = Tem_distance';
    [~, index_dist] = min(Distance); % index_dist(1,:) ���Ƕ�Ӧ�ĵڼ����������ĵ�����ֵ��������� offSubpop(i,:)���壬Ӧ���õ� index_dist(1,:) ���ֲ�����ģ����Ԥɸ
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