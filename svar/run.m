
%% code to generate observed data (dont need to run again since the data used in the paper has been saved)

% nvars = 20;
% X = -0.1*eye(nvars);
% 
% index1 = 1:nvars;
% index2 = 1:nvars;
% 
% pairs = zeros(nvars/2,2);
% theta_true = zeros(nvars+1,1);
% 
% for i = 1:(nvars/2)
% 
%     while(1)
%         r1 = randsample(index1,1);
%         r2 = randsample(index2,1);
%         if (r1 ~= r2)
%             break
%         end
%     end
% 
%     pairs(i,:) = [r1 r2];
% 
%     index1 = index1(index1~=r1);
%     index1 = index1(index1~=r2);
%     index2 = index2(index2~=r1);
%     index2 = index2(index2~=r2);
% 
%     theta_true((i-1)*2+1) = unifrnd(-1,1);
%     theta_true(2*i) = unifrnd(-1,1);
% 
%     X(r1,r2) = theta_true((i-1)*2+1);
%     X(r2,r1) = theta_true(2*i);
% 
% end
% 
% theta_true(end) = 0.1;

%%  SMC ABC with all summaries

% load in saved dataset generated from the model
load('data_GVAR.mat')

extra_args.pairs = pairs;
extra_args.n = length(Y(1,:));
extra_args.nvars = nvars;



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
extra_args.num_params = nvars+1; % this MUST be specified
extra_args.lower = [-1*ones(1,nvars) 0];   %lower lim
extra_args.upper = [1*ones(1,nvars) 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);



smry_func = @(x,extra_args) compute_summaries(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_GVAR(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(Y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
      
save('results_summ.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');   

%%  SMC ABC with all summaries (PILOT)

% load in saved dataset generated from the model
load('data_GVAR.mat')

extra_args.pairs = pairs;
extra_args.n = length(Y(1,:));
extra_args.nvars = nvars;



% the number of particles
N = 1000; 
% target ABC tolerance. (If zero, use acceptance probability stopping criteria)
epsilon_final = 0; 
% tuning parameters for ABC-SMC -- set to good initial defaults
a = 0.5; 
c = 0.01;
% minimum acceptance probability used as stopping criteria if espilon_final = 0 
% if p_acc_min = 0 then use epsilon_final
p_acc_min = 0.10;

% define prior
extra_args.num_params = nvars+1; % this MUST be specified
extra_args.lower = [-1*ones(1,nvars) 0];   %lower lim
extra_args.upper = [1*ones(1,nvars) 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);



smry_func = @(x,extra_args) compute_summaries(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_GVAR(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(Y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
      


save('results_summ_pilot.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');   



%%  SMC ABC with lag1 summary for a particular component in transition matrix

p = 20;  % CHANGE THIS AS NEEDED  (parameter indicator, ranges between 1 and 20)

% load in saved dataset generated from the model
load('data_GVAR.mat')

extra_args.pairs = pairs;
extra_args.n = length(Y(1,:));
extra_args.nvars = nvars;
extra_args.pos = p;



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
extra_args.num_params = nvars+1; % this MUST be specified
extra_args.lower = [-1*ones(1,nvars) 0];   %lower lim
extra_args.upper = [1*ones(1,nvars) 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);



smry_func = @(x,extra_args) compute_summary_lag1(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_GVAR(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(Y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
      

save(['results_summ_p'  num2str(p)  '.mat'],'part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');   


%%  SMC ABC with std summary for std parameter

% load in saved dataset generated from the model
load('data_GVAR.mat')

extra_args.pairs = pairs;
extra_args.n = length(Y(1,:));
extra_args.nvars = nvars;



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
extra_args.num_params = nvars+1; % this MUST be specified
extra_args.lower = [-1*ones(1,nvars) 0];   %lower lim
extra_args.upper = [1*ones(1,nvars) 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);



smry_func = @(x,extra_args) compute_summary_std(x, extra_args); % summary function
sim_func = @(theta,extra_args)simulate_GVAR(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(Y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
      

save('results_summ_std.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');   





%%  SMC ABC continue from pilot with lag1 summary 

p = 20;  % CHANGE THIS AS NEEDED (parameter indicator, ranges between 1 and 20)

% load in saved dataset generated from the model
load('data_GVAR.mat')

extra_args.pairs = pairs;
extra_args.n = length(Y(1,:));
extra_args.nvars = nvars;
extra_args.pos = p;


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
extra_args.num_params = nvars+1; % this MUST be specified
extra_args.lower = [-1*ones(1,nvars) 0];   %lower lim
extra_args.upper = [1*ones(1,nvars) 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);

smry_func_pilot = @(x,extra_args) compute_summaries(x, extra_args); % summary function
smry_func =  @(x,extra_args) compute_summary_lag1(x, extra_args);

sim_func = @(theta,extra_args)simulate_GVAR(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic_continue(Y,part_vals_smc,part_sim_smc(:,p),sim_func,dist_func,dist_func,smry_func,smry_func_pilot,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);

save(['results_summ_continue_p'  num2str(p)  '.mat'],'part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');   



%%  SMC ABC continue from pilot with std summary 

% load in saved dataset generated from the model
load('data_GVAR.mat')

extra_args.pairs = pairs;
extra_args.n = length(Y(1,:));
extra_args.nvars = nvars;


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
extra_args.num_params = nvars+1; % this MUST be specified
extra_args.lower = [-1*ones(1,nvars) 0];   %lower lim
extra_args.upper = [1*ones(1,nvars) 1]; %upper lim

prior_funcs.sampler = @(extra_args) [unifrnd(extra_args.lower,extra_args.upper)];
prior_funcs.trans_f = @(theta,extra_args) [log((theta - extra_args.lower)./(extra_args.upper - theta))]; % logit transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [(extra_args.upper.*exp(theta_trans) + extra_args.lower)./(1 + exp(theta_trans))]; % inverse logit transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);

smry_func_pilot = @(x,extra_args) compute_summaries(x, extra_args); % summary function
smry_func =  @(x,extra_args) compute_summary_std(x, extra_args);

sim_func = @(theta,extra_args)simulate_GVAR(theta,extra_args); % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);

% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic_continue(Y,part_vals_smc,part_sim_smc(:,p),sim_func,dist_func,dist_func,smry_func,smry_func_pilot,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);

save('results_summ_continue_std.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');   





