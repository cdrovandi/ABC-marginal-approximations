function [part_vals, part_sim, part_s, sims,dist_t,p_acc_t,dist_history,sims_history] = smc_abc_rw_generic(y,sim_func,...
                                                 dist_func,smry_func,...
                                                 prior_funcs,extra_args,N,dist_final,...
                                                 a,c,p_acc_min)

%%
% Generic implementation of the SMC ABC replenishment algorithm of Drovandi and Pettitt (2011).  Biometrics.
%%

% compute summary statistics for observations
part_obs = smry_func(y,extra_args);

% values for adaptive steps
num_drop = floor(N*a); %(usually a = 0.5) so keep half
num_keep = N-num_drop;
mcmc_trials = 5;

% Initialise particle data structures
part_vals = zeros(N,extra_args.num_params);
part_s = zeros(N,1);
part_sim = zeros(N, length(part_obs));

% initial prior rejection algorithm
parfor i = 1:N
    % sample prior 
    part_vals(i,:) = prior_funcs.sampler(extra_args);
    % simulate model
    y_s = sim_func(part_vals(i,:), extra_args);
    part_sim(i,:) = smry_func(y_s, extra_args);
    % evaluate the discrepancy metric    
    part_s(i) = dist_func(part_obs,part_sim(i,:), extra_args);
end

sims = N;

dist_history = max(part_s);
sims_history = N;

% transform the parameters
for i=1:N
    part_vals(i,:) = prior_funcs.trans_f(part_vals(i,:),extra_args);
end

% sort the particles
[part_s,ix] = sort(part_s); 
part_vals = part_vals(ix,:); 
part_sim = part_sim(ix,:);

% determine next disprepacy threshold
dist_max = part_s(N);
dist_next = part_s(num_keep);
dist_max,dist_next,dist_final
dist_t = dist_next;
p_acc_t = 0;

% interate toward target discrepancy
while (dist_max > dist_final)
    fprintf('****************************** \n');
    % compute the covariance matrix (of particles that remain) required 
    % for the Independent MH move step
    cov_matrix = (2.38^2)*cov(part_vals(1:num_keep,:))/size(part_vals,2);
    
    % resample
    r = randsample(num_keep, N-num_keep, 'true');
    part_vals((num_keep+1):N, :) = part_vals(r,:);
    part_s((num_keep+1):N) = part_s(r);
    part_sim((num_keep+1):N, :) = part_sim(r,:);
    
    i_acc = zeros(N-num_keep,1);
    sims_mcmc = zeros(N-num_keep,1);
    
    % trial MCMC steps 
    parfor i = (num_keep+1):N
        for r = 1:mcmc_trials
            % Gaussian random walk
            part_vals_prop = mvnrnd(part_vals(i,:),cov_matrix);
            prior_curr = prior_funcs.pdf(part_vals(i,:),extra_args);
            prior_prop = prior_funcs.pdf(part_vals_prop,extra_args);
            % early rejection (assumes symmetric proposals)
            if (isnan(prior_prop/prior_curr) || rand > prior_prop/prior_curr)  
                continue;
            end

            prop = prior_funcs.trans_finv(part_vals_prop,extra_args);
            
            y_s = sim_func(prop, extra_args);
            part_sim_prop = smry_func(y_s,extra_args);
            dist_prop = dist_func(part_obs,part_sim_prop,extra_args);
            
            sims_mcmc(i-num_keep) = sims_mcmc(i-num_keep)+1;
            % ABC part of the acceptance probability
            
            %Accept a particle if it is within the target distance.
            if (dist_prop <= dist_next) 
                % then the metropolis-hastings ratio is satisfied
                part_vals(i,:) = part_vals_prop; 
                part_s(i) = dist_prop; 
                part_sim(i,:) = part_sim_prop;
                i_acc(i-num_keep) = i_acc(i-num_keep) + 1;
            end
        end
    end
    
    % determine number of MCMC iterations to perfrom 
    acc_rate = sum(i_acc)/(mcmc_trials*(N-num_keep));
    mcmc_iters =   floor(log(c)/log(1-acc_rate)+1);
    fprintf('Total number of mcmc moves for current target is %d, number remaining is %d\n',mcmc_iters,mcmc_iters-mcmc_trials);
    
    % move step
    parfor i = (num_keep+1):N
        for r = 1:(mcmc_iters-mcmc_trials)
            % Gaussian random walk
            part_vals_prop = mvnrnd(part_vals(i,:),cov_matrix);
            prior_curr = prior_funcs.pdf(part_vals(i,:), extra_args);
            prior_prop = prior_funcs.pdf(part_vals_prop, extra_args);
            % early rejection (assumes symmetric proposal)
            if (isnan(prior_prop/prior_curr) || rand > prior_prop/prior_curr) 
                continue;
            end
            
            prop = prior_funcs.trans_finv(part_vals_prop,extra_args);
            y_s = sim_func(prop, extra_args);
            part_sim_prop = smry_func(y_s,extra_args);
            dist_prop = dist_func(part_obs,part_sim_prop,extra_args);
            
            sims_mcmc(i-num_keep) = sims_mcmc(i-num_keep)+1;
            
            % ABC part of the acceptance probability
            if (dist_prop <= dist_next) 
                % then the metropolis-hastings ratio is satisfied
                part_vals(i,:) = part_vals_prop; 
                part_s(i) = dist_prop; 
                part_sim(i,:) = part_sim_prop;
                i_acc(i-num_keep) = i_acc(i-num_keep) + 1;
            end
        end
    end
    
    num_mcmc_iters = max(0, mcmc_iters - mcmc_trials) + mcmc_trials;
    p_acc = sum(i_acc)/(num_mcmc_iters*(N-num_keep));
    
    fprintf('MCMC acceptance probability was %f\n',p_acc);
    
    sims = sims + sum(sims_mcmc);
    mcmc_trials = ceil(mcmc_iters/2);
    
    % compute number of unique particles
    fprintf('The number of unique particles is %d\n',length(unique(part_vals(:,1))));
    
    % compute the next distance and maximum distance
    % sort the particles
    [part_s,ix] = sort(part_s); part_vals = part_vals(ix,:); part_sim = part_sim(ix,:);
    
    % if most of the particles are under the final target then don't
    % drop them
    if (sum((part_s > dist_final)) < num_drop)
        num_drop = sum((part_s > dist_final));
        num_keep = N-num_drop;
    end
   
    % to return information about convergence
    dist_t = dist_next;
    p_acc_t = p_acc;

    dist_max = part_s(N);
    dist_next = part_s(num_keep);
    % check to see if we reach desired tolerance at next iteration
    if (dist_next < dist_final)
        dist_next = dist_final;
    end
    
    fprintf('The next distance is %f and the maximum distance is %f and the number to drop is %d\n',dist_next,dist_max,num_drop);
    fprintf('The number of sims is %d\n',sims);
    
    dist_history = [dist_history max(part_s)];
    sims_history = [sims_history sims];
   
    if (p_acc < p_acc_min)
        fprintf('Getting out as MCMC acceptance rate is below acceptable threshold\n');
        break;
    end

end
for i=1:N
    part_vals(i,:) = prior_funcs.trans_finv(part_vals(i,:),extra_args);
end
            
