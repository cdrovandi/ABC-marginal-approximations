function f = summary_post_mean_variance(x,extra_args)

	% summary statistics based off laplace approximation

    comp = extra_args.comp;
    [phat(1) phat(2)] = normfit(x);
    
    options = optimoptions(@fminunc,'Display','off');
	
	[post_mode,~,~,~,~,hessian] = fminunc(@(theta_trans) normlike([theta_trans(1) sqrt(exp(theta_trans(2)))], x) - log(normpdf(theta_trans(1), extra_args.mu0, sqrt(extra_args.phi0))) - (-extra_args.a - 1)*theta_trans(2) + extra_args.b/exp(theta_trans(2)) - theta_trans(2), [phat(1) log(phat(2)^2)], options);
    
    post_std = diag(sqrt(inv(hessian)));
    
    f = [post_mode(comp) post_std(comp)];
    
end