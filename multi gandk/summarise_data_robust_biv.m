function [summ] = summarise_data_robust_biv(data, extra_args)
%%
% Robust summary statistics for bivariate g-and-k distribution
%%
    num_obs = extra_args.n;

    % location, scale, skewness and kurtosis
    summ(1:4) = summarise_data_robust(data(:,1));
    summ(5:8) = summarise_data_robust(data(:,2));
    
    r1 = tiedrank(data(:,1));
    r2 = tiedrank(data(:,2));
    z1 = norminv(r1./(num_obs+1));
    z2 = norminv(r2./(num_obs+1));
    
    summ(9) = sum(z1.*z2)/sum(z1.^2);
    
end