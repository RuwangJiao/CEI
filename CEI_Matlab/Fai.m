function Fai=Fai(z)
    global integral_data
    global Dmodel
    global index_dist
    Num      = integral_data.Num;
    x        = integral_data.x;
    PI_total = 1;
    for i = 1:Num
        [offspring_surrogate, dy, mse] = predictor(x, Dmodel(index_dist, i+1).dmodel);
        y_hat                          = offspring_surrogate;
        SSqr                           = mse;
        probImp1                       = normcdf((z-y_hat)./sqrt(SSqr));
        PI_total                       = PI_total.*probImp1;
    end
    Fai = PI_total;
end