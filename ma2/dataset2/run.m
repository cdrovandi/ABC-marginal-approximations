
                              
           


%% set-up ABC-SMC 
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
extra_args.num_params = 2; % this MUST be specified

prior_funcs.sampler = @(extra_args) prior_sampler(extra_args);
prior_funcs.trans_f = @(theta,extra_args) [theta]; % no transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans]; % no transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prior_pdf(theta_trans,extra_args);


load('data.mat')
extra_args.T = length(y);
extra_args.K = 3;
extra_args.reps = 1;

smry_func = @(x,extra_args) summStat_ma2_acf(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_ma2(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
       



%% set-up ABC-SMC - pilot
% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.3;



% define prior
extra_args.num_params = 2; % this MUST be specified

prior_funcs.sampler = @(extra_args) prior_sampler(extra_args);
prior_funcs.trans_f = @(theta,extra_args) [theta]; % no transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans]; % no transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prior_pdf(theta_trans,extra_args);


load('data.mat')
extra_args.T = length(y);
extra_args.K = 3;
extra_args.reps = 1;

smry_func = @(x,extra_args) summStat_ma2_acf(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_ma2(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_pilot.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
                                   


%% set-up ABC-SMC - with variance and first autocovariance
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
extra_args.num_params = 2; % this MUST be specified

prior_funcs.sampler = @(extra_args) prior_sampler(extra_args);
prior_funcs.trans_f = @(theta,extra_args) [theta]; % no transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans]; % no transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prior_pdf(theta_trans,extra_args);


load('data.mat')
extra_args.T = length(y);
extra_args.reps = 1;

smry_func = @(x,extra_args) summStat_ma2_acf1(x, extra_args); % summary function - this time just use the variance and first autocovariance
sim_func = @(theta,extra_args)simulate_ma2(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_autocov1.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
       


%% set-up ABC-SMC - with second autocovariance
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
extra_args.num_params = 2; % this MUST be specified

prior_funcs.sampler = @(extra_args) prior_sampler(extra_args);
prior_funcs.trans_f = @(theta,extra_args) [theta]; % no transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans]; % no transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prior_pdf(theta_trans,extra_args);


load('data.mat')
extra_args.T = length(y);
extra_args.reps = 1;

smry_func = @(x,extra_args) summStat_ma2_acf2(x, extra_args); % summary function - this time just second autocovariance
sim_func = @(theta,extra_args)simulate_ma2(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_autocov2.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
       




%% set-up ABC-SMC - continue from pilot with variance and first autovariance
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
extra_args.num_params = 2; % this MUST be specified

prior_funcs.sampler = @(extra_args) prior_sampler(extra_args);
prior_funcs.trans_f = @(theta,extra_args) [theta]; % no transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans]; % no transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prior_pdf(theta_trans,extra_args);


load('data.mat')
extra_args.T = length(y);
extra_args.K = 3; 
extra_args.reps = 1;

smry_func_pilot = @(x,extra_args) summStat_ma2_acf(x, extra_args); % summary function
smry_func = @(x,extra_args) summStat_ma2_acf1(x, extra_args);

sim_func = @(theta,extra_args)simulate_ma2(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic_continue(y,part_vals_smc,sim_func,dist_func,dist_func,smry_func,smry_func_pilot,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_autocov1_continue.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
        


%% set-up ABC-SMC - continue from pilot with second autovariance
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
extra_args.num_params = 2; % this MUST be specified

prior_funcs.sampler = @(extra_args) prior_sampler(extra_args);
prior_funcs.trans_f = @(theta,extra_args) [theta]; % no transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans]; % no transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prior_pdf(theta_trans,extra_args);


load('data.mat')
extra_args.T = length(y);
extra_args.K = 3; 
extra_args.reps = 1;

smry_func_pilot = @(x,extra_args) summStat_ma2_acf(x, extra_args); % summary function
smry_func = @(x,extra_args) summStat_ma2_acf2(x, extra_args);

sim_func = @(theta,extra_args)simulate_ma2(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic_continue(y,part_vals_smc,sim_func,dist_func,dist_func,smry_func,smry_func_pilot,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_autocov2_continue.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
        
