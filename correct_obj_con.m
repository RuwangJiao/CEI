function [h, num_inequal, num_equal, fit] = correct_obj_con(h1, g, fit)
    h                           =h1(:,1:end/2);
    num_inequal                 = size(g,2);
    num_equal                   = size(h,2);
    fit(:, end-num_equal+1:end) = [];
end