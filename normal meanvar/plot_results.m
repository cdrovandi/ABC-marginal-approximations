
load('results_summ.mat');
part_vals = part_vals_smc;
load('results_summ1.mat');
part_vals1 = part_vals_smc;
load('results_summ2.mat');
part_vals2 = part_vals_smc;


load('results_summ_mean1.mat')
part_vals_mean1 = part_vals_smc;
load('results_summ_mean2.mat')
part_vals_mean2 = part_vals_smc;


figure;
[f,xi] = ksdensity(part_vals(:,1));
plot(xi,f,'k','LineWidth',2);
hold on
[f,xi] = ksdensity(part_vals_mean1(:,1));
plot(xi,f,'r-.','LineWidth',2)
[f,xi] = ksdensity(part_vals1(:,1));
plot(xi,f,'b--','LineWidth',2)
xlim([-1 1]);
xlabel('\mu','FontSize',16);
ylabel('density','FontSize',16);
legend('sufficient', 'first moment', 'first two moments')

figure;
[f,xi] = ksdensity(part_vals(:,2));
plot(xi,f,'k','LineWidth',2);
hold on
[f,xi] = ksdensity(part_vals_mean2(:,2));
plot(xi,f,'r-.','LineWidth',2)
[f,xi] = ksdensity(part_vals2(:,2));
plot(xi,f,'b--','LineWidth',2)
xlim([0.4 1.6]);
xlabel('\phi','FontSize',16);
ylabel('density','FontSize',16);
legend('sufficient', 'first moment', 'first two moments')



