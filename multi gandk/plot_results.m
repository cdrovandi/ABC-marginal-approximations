

load('results_summ.mat')
part_vals_all = part_vals_smc;
part_summ_all = part_sim_smc;

load('results_summ_pilot.mat')
part_vals_pilot = part_vals_smc;
part_summ_pilot = part_sim_smc;

load('results_summ_corr.mat')
part_vals_corr = part_vals_smc;
part_summ_corr = part_sim_smc;

load('results_summ_corr_continue.mat')
part_vals_corr_cont = part_vals_smc;
part_summ_corr_cont = part_sim_smc;


figure;
[f,xi] = ksdensity(part_vals_all(:,9));
plot(xi,f,'k','LineWidth',2);
hold on;
[f,xi] = ksdensity(part_vals_corr(:,9));
plot(xi,f,'b--','LineWidth',2);
plot(0.6,0,'kx','MarkerSize',12)
xlabel('r','FontSize',16);
ylabel('density','FontSize',16);
legend('\pi(r|S)', '\pi(r|S_r)');


figure;
[f,xi] = ksdensity(part_summ_all(:,9));
plot(xi,f,'k','LineWidth',2);
hold on;
[f,xi] = ksdensity(part_summ_pilot(:,9));
plot(xi,f,'b--','LineWidth',2);
[f,xi] = ksdensity(part_summ_corr);
plot(xi,f,'g-.','LineWidth',2);
[f,xi] = ksdensity(part_summ_corr_cont);
plot(xi,f,'r.','LineWidth',2);
xlabel('S','FontSize',16);
ylabel('density','FontSize',16);
legend('all summs', 'pilot (all summs)', 'S', 'pilot + S')

