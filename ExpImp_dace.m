function ExpImp1 = ExpImp_dace(x, y_min, dmodel)
    [offspring_surrogate, dy, mse] = predictor(x, dmodel);
    y_hat                          = offspring_surrogate;
    SSqr                           = mse;
    %  Expected   improvement
    if SSqr == 0
        ExpImp1 = 0;  
    else
    % ei_termone = (y_min-y_hat)*(0.5 + 0.5*erf((1/sqrt(2))*((y_min  - y_hat)/sqrt(abs(SSqr)))));
    % ei_termtwo = sqrt(abs(SSqr))*(1/sqrt(2*pi))*exp( - (1/2)*((y_min - y_hat)^2/SSqr));
    % ExpImp1 = ei_termone + ei_termtwo;
        ei_termone = (y_min - y_hat)*normcdf((y_min  - y_hat)/sqrt(abs(SSqr)));
        ei_termtwo = sqrt(abs(SSqr))*normpdf((y_min  - y_hat)/sqrt(abs(SSqr)));
        ExpImp1    = ei_termone + ei_termtwo;
    end
    % ExpImp1=-ExpImp1;
end
