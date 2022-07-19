function f = summStat_mvn_component(x, extra_args)
    Sigma_P = extra_args.Sigma_P;
    Sigma = extra_args.Sigma;
    
    mu = Sigma_P*inv(Sigma)*x;

    f = mu(extra_args.p);
    
end