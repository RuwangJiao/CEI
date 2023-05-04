function [h, num_inequal, num_equal, fit] = correct_obj_con_equ(h1, g, fit)
    h           = h1;
    num_inequal = size(g, 2);
    num_equal   = floor(size(h, 2)/2);
end