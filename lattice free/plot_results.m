
load('data_simulated/data.mat','y','init_data');
smry_func = @(x,extra_args) summStat_pair_corr(x,extra_args);
extra_args.PC_dr = 50;
obs_summ = summStat_num_cells(y,extra_args);

load('results_summ.mat')
part_vals_all = part_vals_smc;
part_summ_all = part_sim_smc;

load('results_summ_pilot.mat')
part_vals_pilot = part_vals_smc;
part_summ_pilot = part_sim_smc;

load('results_summ_numcells.mat')
part_vals_numcells = part_vals_smc;
part_summ_numcells = part_sim_smc;

load('results_summ_numcells_continue.mat')
part_vals_numcells_cont = part_vals_smc;
part_summ_numcells_cont = part_sim_smc;


figure;
[f,xi] = ksdensity(part_vals_all(:,2));
plot(xi,f,'k','LineWidth',2);
hold on;
[f,xi] = ksdensity(part_vals_pilot(:,2));
plot(xi,f,'b--','LineWidth',2);
[f,xi] = ksdensity(part_vals_numcells(:,2));
plot(xi,f,'g-.','LineWidth',2);
[f,xi] = ksdensity(part_vals_numcells_cont(:,2));
plot(xi,f,'r.','LineWidth',2);
plot(0.04,0,'kx','MarkerSize',12)
xlim([0.028 0.054]);
xlabel('\rho','FontSize',16);
ylabel('density','FontSize',16);
legend('S', 'pilot S', 'S_p', 'pilot + S_p')


figure;
[f,xi] = ksdensity(part_summ_all(:,5));
plot(xi,f,'k','LineWidth',2);
hold on;
[f,xi] = ksdensity(part_summ_pilot(:,5));
plot(xi,f,'b--','LineWidth',2);
[f,xi] = ksdensity(part_summ_numcells);
plot(xi,f,'g-.','LineWidth',2);
[f,xi] = ksdensity(part_summ_numcells_cont);
plot(xi,f,'r.','LineWidth',2);
plot(obs_summ,0,'kx','MarkerSize',12)
xlim([1270 1320]);
xlabel('S_p','FontSize',16);
ylabel('density','FontSize',16);
legend('S', 'pilot S', 'S_p', 'pilot + S_p')

