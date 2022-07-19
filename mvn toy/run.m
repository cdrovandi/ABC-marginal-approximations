
%% code to simulate data

% covariance matrix of data
d = 10;
Sigma = zeros(d,d);
for i = 1:d
    for j = 1:d
        if (i == j)
            Sigma(i,i) = 1;
        else
            Sigma(i,j) = 0.9^(abs(i-j));
            Sigma(j,i) = Sigma(i,j);
        end
        
    end
end

y = mvnrnd(zeros(d,1), Sigma)';

% posterior mean and covariance
Sigma_P = inv(eye(d) + inv(Sigma));
mu_P = Sigma_P*inv(Sigma)*y;


%% running of SMC ABC algorithm

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

p = 1; % component mean we are trying to estimate

% define prior
extra_args.num_params = 10; % this MUST be specified

prior_funcs.sampler = @(extra_args) [normrnd(0,1,1,extra_args.num_params)];
prior_funcs.trans_f = @(theta,extra_args) [theta]; % no transform
prior_funcs.trans_finv = @(theta_trans,extra_args) [theta_trans]; % no transform
% the prior pdf must be the density for the transformed space (ie the space that MCMC samples over in the move step)
prior_funcs.pdf = @(theta_trans,extra_args) [prod(normpdf(theta_trans,0,1))];


load('data.mat')
extra_args.Sigma = Sigma;
extra_args.Sigma_P = Sigma_P;
extra_args.p = p;

smry_func = @(x,extra_args) summStat_mvn_component(x, extra_args); % summary function
sim_func = @(theta,extra_args)[mvnrnd(theta, extra_args.Sigma)']; % simulation function

%discrepancy metric
dist_func = @(sim_smry,obs_smry,extra_args) sum((sim_smry-obs_smry).^2);


% run SMC sampler
[part_vals_smc, part_sim_smc, part_s_smc, smc_sims, smc_epsilon_t, smc_p_acc_t] = smc_abc_rw_generic(y,sim_func,dist_func,smry_func,prior_funcs,extra_args,N,epsilon_final,a,c,p_acc_min);
       
save('results_summ1.mat','part_vals_smc','part_sim_smc','part_s_smc','smc_sims','smc_epsilon_t','smc_p_acc_t');                    
               




