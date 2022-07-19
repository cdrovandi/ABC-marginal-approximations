function x = simulate_biv_gk(parms,extra_args)
%%
% method to simulate from the bivariate g and k distribution
% a1 = parms(1); b1 = parms(2); g1 = parms(3); k1 = parms(4);
% a2 = parms(5); b2 = parms(6); g2 = parms(7); k2 = parms(8);
% rho = parms(9);
% n - is the number of observations
%%

    n = extra_args.n;

    rho = parms(9);
    z1 = randn(1,n);
    z2 = rho*z1 + sqrt(1-rho^2)*randn(1,n);
    
    x(:,1) = fun_gk(parms(1:4),z1);
    x(:,2) = fun_gk(parms(5:8),z2);
end