function [summ] = summarise_data_corr(data, extra_args)
%%
% Robust summary statistics for bivariate g-and-k distribution (just the rank correlation)
%%
    num_obs = extra_args.n;
    
    r1 = tiedrank(data(:,1));
    r2 = tiedrank(data(:,2));
    z1 = norminv(r1./(num_obs+1));
    z2 = norminv(r2./(num_obs+1));
    
    summ = sum(z1.*z2)/sum(z1.^2);
    
end