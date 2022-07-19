





%% SMC ABC bivariate all summaries 

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

% define prior
extra_args.num_params = 9; % this MUST be specified
extra_args.lower = [0 0 0 0 0 0 0 0 -1];   %lower lim
extra_args.upper = [10 10 10 10 10 10 10 10 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);


% Simulated Data - Bivariate Dataset correlation is 0.6 - moderate
y = load('bivariate simulated dataset corr 06.txt');


extra_args.n = length(y(:,1));

smry_func = @(x,extra_args) summarise_data_robust_biv(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_biv_gk(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   


%% SMC ABC bivariate all summaries - PILOT

% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.12;

% define prior
extra_args.num_params = 9; % this MUST be specified
extra_args.lower = [0 0 0 0 0 0 0 0 -1];   %lower lim
extra_args.upper = [10 10 10 10 10 10 10 10 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);


% Simulated Data - Bivariate Dataset correlation is 0.6 - moderate
y = load('bivariate simulated dataset corr 06.txt');


extra_args.n = length(y(:,1));

smry_func = @(x,extra_args) summarise_data_robust_biv(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_biv_gk(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_pilot.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   


%% SMC ABC bivariate just correlation summary 

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

% define prior
extra_args.num_params = 9; % this MUST be specified
extra_args.lower = [0 0 0 0 0 0 0 0 -1];   %lower lim
extra_args.upper = [10 10 10 10 10 10 10 10 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);


% Simulated Data - Bivariate Dataset correlation is 0.6 - moderate
y = load('bivariate simulated dataset corr 06.txt');


extra_args.n = length(y(:,1));

smry_func = @(x,extra_args) summarise_data_corr(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_biv_gk(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_corr.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
       
	  
%% SMC ABC bivariate continue from pilot with correlation summary 
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

load('results_summ_pilot.mat');

extra_args.dist_pilot = max(part_s_smc);

% define prior
extra_args.num_params = 9; % this MUST be specified
extra_args.lower = [0 0 0 0 0 0 0 0 -1];   %lower lim
extra_args.upper = [10 10 10 10 10 10 10 10 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)]; 
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);

% Simulated Data - Bivariate Dataset correlation is 0.6 - moderate
y = load('bivariate simulated dataset corr 06.txt');
extra_args.n = length(y(:,1));

smry_func_pilot = @(x,extra_args) summarise_data_robust_biv(x, extra_args);
smry_func = @(x,extra_args) summarise_data_corr(x, extra_args);
sim_func = @(theta,extra_args)simulate_biv_gk(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum(((sim_smry-obs_smry)./obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic_continue(y,part_vals_smc,part_sim_smc(:,9),sim_func,dist_func,dist_func,smry_func,smry_func_pilot,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_corr_continue.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
        
	   
	   