
                              
           


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
extra_args.num_params = 3; % this MUST be specified
extra_args.lower = [0 0 0];   %lower lim
extra_args.upper = [10 0.1 20]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)]; 
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);


load('data_simulated/data.mat','y','init_data');
smry_func = @(x,extra_args) summStat_pair_corr(x,extra_args);
extra_args.PC_dr = 50;

sim_func = @(theta,extra_args)simulate_lattice_free_cell(theta,extra_args);
extra_args.init_data = init_data;

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum(((sim_smry-obs_smry)./obs_smry).^2);


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
p_acc_min = 0.15;



% define prior
extra_args.num_params = 3; % this MUST be specified
extra_args.lower = [0 0 0];   %lower lim
extra_args.upper = [10 0.1 20]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)]; 
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);


load('data_simulated/data.mat','y','init_data');
smry_func = @(x,extra_args) summStat_pair_corr(x,extra_args);
extra_args.PC_dr = 50;

sim_func = @(theta,extra_args)simulate_lattice_free_cell(theta,extra_args);
extra_args.init_data = init_data;

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum(((sim_smry-obs_smry)./obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_pilot.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
                                   

             



                        
%% set-up ABC-SMC - number of cells summary statisticss
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
extra_args.num_params = 3; % this MUST be specified
extra_args.lower = [0 0 0];   %lower lim
extra_args.upper = [10 0.1 20]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)]; 
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);


load('data_simulated/data.mat','y','init_data');
smry_func = @(x,extra_args) summStat_num_cells(x,extra_args);
extra_args.PC_dr = 50;

sim_func = @(theta,extra_args)simulate_lattice_free_cell(theta,extra_args);
extra_args.init_data = init_data;

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum(((sim_smry-obs_smry)./obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_numcells.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
                                   



%% set-up ABC-SMC - continue
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
extra_args.num_params = 3; % this MUST be specified
extra_args.lower = [0 0 0];   %lower lim
extra_args.upper = [10 0.1 20]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)]; 
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);


load('data_simulated/data.mat','y','init_data');
smry_func_pilot = @(x,extra_args) summStat_pair_corr(x,extra_args);
smry_func = @(x,extra_args) summStat_num_cells(x,extra_args);
extra_args.PC_dr = 50;

sim_func = @(theta,extra_args)simulate_lattice_free_cell(theta,extra_args);
extra_args.init_data = init_data;

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum(((sim_smry-obs_smry)./obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic_continue(y,part_vals_smc,part_sim_smc(:,5),sim_func,dist_func,dist_func,smry_func,smry_func_pilot,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ_numcells_continue.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
                                   
        