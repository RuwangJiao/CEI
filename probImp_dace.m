function probImp1 = probImp_dace(x, improv_y, dmodel)
    [offspring_surrogate, dy, mse] = predictor(x, dmodel);
    y_hat                          = offspring_surrogate;
    SSqr                           = mse;
    probImp1                       = normcdf((improv_y - y_hat)/sqrt(abs(SSqr)));
end