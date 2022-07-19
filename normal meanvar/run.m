

n = 100;
mu = 0; phi = 1;
y = normrnd(mu,sqrt(phi),n,1);




%% SMC ABC all data 

load('data.mat');
extra_args.num_params = 2; % this MUST be specified
extra_args.n = n;
extra_args.mu0 = 0;
extra_args.phi0 = 1;
extra_args.a = 1;
extra_args.b = 1;


% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.01;

prior_funcs.sampler = @(extra_args) [normrnd(extra_args.mu0, sqrt(extra_args.phi0)) 1/gamrnd(extra_args.a, 1/extra_args.b)];
prior_funcs.trans_f = @(theta,extra_args) [theta(1) log(theta(2))]; 
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans(1) exp(theta_trans(2))]; 
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args)  [normpdf(theta_trans(1), extra_args.mu0, sqrt(extra_args.phi0))*exp(theta_trans(2))^(-extra_args.a - 1)*exp(-extra_args.b/exp(theta_trans(2)))*exp(theta_trans(2))];

smry_func = @(x,extra_args) [mean(x) std(x)]; % summary function
sim_func = @(theta,extra_args) normrnd(theta(1), sqrt(theta(2)), extra_args.n, 1); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       

save('results_summ.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   


%% SMC ABC - just posterior mean and variance of the first parameter as summaries


load('data.mat');
extra_args.num_params = 2; % this MUST be specified
extra_args.n = n;
extra_args.mu0 = 0;
extra_args.phi0 = 1;
extra_args.a = 1;
extra_args.b = 1;
extra_args.comp = 1;


% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.01;

prior_funcs.sampler = @(extra_args) [normrnd(extra_args.mu0, sqrt(extra_args.phi0)) 1/gamrnd(extra_args.a, 1/extra_args.b)];
prior_funcs.trans_f = @(theta,extra_args) [theta(1) log(theta(2))]; 
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans(1) exp(theta_trans(2))]; 
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args)  [normpdf(theta_trans(1), extra_args.mu0, sqrt(extra_args.phi0))*exp(theta_trans(2))^(-extra_args.a - 1)*exp(-extra_args.b/exp(theta_trans(2)))*exp(theta_trans(2))];

smry_func = @(x,extra_args) summary_post_mean_variance(x,extra_args); % summary function
sim_func = @(theta,extra_args) normrnd(theta(1), sqrt(theta(2)), extra_args.n, 1); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       

save('results_summ1.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
              


%% SMC ABC - just posterior mean of the first parameter as summaries


load('data.mat');
extra_args.num_params = 2; % this MUST be specified
extra_args.n = n;
extra_args.mu0 = 0;
extra_args.phi0 = 1;
extra_args.a = 1;
extra_args.b = 1;
extra_args.comp = 1;


% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.01;

prior_funcs.sampler = @(extra_args) [normrnd(extra_args.mu0, sqrt(extra_args.phi0)) 1/gamrnd(extra_args.a, 1/extra_args.b)];
prior_funcs.trans_f = @(theta,extra_args) [theta(1) log(theta(2))]; 
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans(1) exp(theta_trans(2))]; 
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args)  [normpdf(theta_trans(1), extra_args.mu0, sqrt(extra_args.phi0))*exp(theta_trans(2))^(-extra_args.a - 1)*exp(-extra_args.b/exp(theta_trans(2)))*exp(theta_trans(2))];

smry_func = @(x,extra_args) summary_post_mean_variance(x,extra_args); % summary function
sim_func = @(theta,extra_args) normrnd(theta(1), sqrt(theta(2)), extra_args.n, 1); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry(1)-obs_smry(1)).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       

save('results_summ_mean1.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
              


%% SMC ABC - just posterior mean and variance of the second parameter as summaries


load('data.mat');
extra_args.num_params = 2; % this MUST be specified
extra_args.n = n;
extra_args.mu0 = 0;
extra_args.phi0 = 1;
extra_args.a = 1;
extra_args.b = 1;
extra_args.comp = 2;


% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.01;

prior_funcs.sampler = @(extra_args) [normrnd(extra_args.mu0, sqrt(extra_args.phi0)) 1/gamrnd(extra_args.a, 1/extra_args.b)];
prior_funcs.trans_f = @(theta,extra_args) [theta(1) log(theta(2))]; 
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans(1) exp(theta_trans(2))]; 
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args)  [normpdf(theta_trans(1), extra_args.mu0, sqrt(extra_args.phi0))*exp(theta_trans(2))^(-extra_args.a - 1)*exp(-extra_args.b/exp(theta_trans(2)))*exp(theta_trans(2))];

smry_func = @(x,extra_args) summary_post_mean_variance(x,extra_args); % summary function
sim_func = @(theta,extra_args) normrnd(theta(1), sqrt(theta(2)), extra_args.n, 1); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       

save('results_summ2.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
              


%% SMC ABC - just posterior mean of the second parameter as summaries


load('data.mat');
extra_args.num_params = 2; % this MUST be specified
extra_args.n = n;
extra_args.mu0 = 0;
extra_args.phi0 = 1;
extra_args.a = 1;
extra_args.b = 1;
extra_args.comp = 2;


% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.01;

prior_funcs.sampler = @(extra_args) [normrnd(extra_args.mu0, sqrt(extra_args.phi0)) 1/gamrnd(extra_args.a, 1/extra_args.b)];
prior_funcs.trans_f = @(theta,extra_args) [theta(1) log(theta(2))]; 
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans(1) exp(theta_trans(2))]; 
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args)  [normpdf(theta_trans(1), extra_args.mu0, sqrt(extra_args.phi0))*exp(theta_trans(2))^(-extra_args.a - 1)*exp(-extra_args.b/exp(theta_trans(2)))*exp(theta_trans(2))];

smry_func = @(x,extra_args) summary_post_mean_variance(x,extra_args); % summary function
sim_func = @(theta,extra_args) normrnd(theta(1), sqrt(theta(2)), extra_args.n, 1); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry(1)-obs_smry(1)).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       

save('results_summ_mean2.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
              


