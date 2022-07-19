

load('results_summ.mat')
part_vals_all = part_vals_smc;
part_summ_all = part_sim_smc;

load('results_summ_pilot.mat')
part_vals_pilot = part_vals_smc;
part_summ_pilot = part_sim_smc;

load('results_summ_autocov1.mat')
part_vals_autocov1 = part_vals_smc;
part_summ_autocov1 = part_sim_smc;

load('results_summ_autocov1_continue.mat')
part_vals_autocov1_cont = part_vals_smc;
part_summ_autocov1_cont = part_sim_smc;



figure;
subaxis('MT',0.02,'MB',0.05,'ML',0.05,'MR',0.02,'PL',0.02,'PR',0.01,'Pt',0.015,'PB',0.1);

subaxis(1,2,1);

[f,xi] = ksdensity(part_vals_all(:,1));
plot(xi,f,'k','LineWidth',2);
hold on;
[f,xi] = ksdensity(part_vals_pilot(:,1));
plot(xi,f,'b--','LineWidth',2);
[f,xi] = ksdensity(part_vals_autocov1(:,1));
plot(xi,f,'g-.','LineWidth',2);
[f,xi] = ksdensity(part_vals_autocov1_cont(:,1));
plot(xi,f,'r.','LineWidth',2);
plot(0.9,0,'kx','MarkerSize',12)
xlim([0.38 1]);
xlabel('\theta_1','FontSize',16);
ylabel('density','FontSize',16);
legend('S_1S_2S_3', 'pilot', 'S_1S_2', 'pilot + S_1S_2')




load('results_summ.mat')
part_vals_all = part_vals_smc;
part_summ_all = part_sim_smc;

load('results_summ_pilot.mat')
part_vals_pilot = part_vals_smc;
part_summ_pilot = part_sim_smc;

load('results_summ_autocov2.mat')
part_vals_autocov2 = part_vals_smc;
part_summ_autocov2 = part_sim_smc;

load('results_summ_autocov2_continue.mat')
part_vals_autocov2_cont = part_vals_smc;
part_summ_autocov2_cont = part_sim_smc;



subaxis(1,2,2);

[f,xi] = ksdensity(part_vals_all(:,2));
plot(xi,f,'k','LineWidth',2);
hold on;
[f,xi] = ksdensity(part_vals_pilot(:,2));
plot(xi,f,'b--','LineWidth',2);
[f,xi] = ksdensity(part_vals_autocov2(:,2));
plot(xi,f,'g-.','LineWidth',2);
[f,xi] = ksdensity(part_vals_autocov2_cont(:,2));
plot(xi,f,'r.','LineWidth',2);
plot(-0.05,0,'kx','MarkerSize',12)
xlim([-0.17 0.07]);
xlabel('\theta_2','FontSize',16);
ylabel('density','FontSize',16);
legend('S_1S_2S_3', 'pilot', 'S_3', 'pilot + S_3')







