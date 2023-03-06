
% produces figure in paper based on the generated results

load('data_GVAR.mat')

load('results_summ.mat')
part_vals_summ = part_vals_smc;

load('results_summ_pilot.mat')
part_vals_summ_pilot = part_vals_smc;


t = tiledlayout('flow','TileSpacing','compact','Padding','tight');

for p = 1:20
    nexttile;

    load(['results_summ_p'  num2str(p)  '.mat']);
    part_vals_single = part_vals_smc;

    load(['results_summ_continue_p'  num2str(p)  '.mat']);
    part_vals_single_continue = part_vals_smc;


    [f,xi] = ksdensity(part_vals_summ(:,p));
    plot(xi,f,'k','LineWidth',2);
    hold on;
    [f,xi] = ksdensity(part_vals_summ_pilot(:,p));
    plot(xi,f,'b--','LineWidth',2);
    [f,xi] = ksdensity(part_vals_single(:,p));
    plot(xi,f,'g-.','LineWidth',2);
    [f,xi] = ksdensity(part_vals_single_continue(:,p));
    plot(xi,f,'r.','LineWidth',2);
    plot(theta_true(p),0,'kx');



end

p=21;
load('results_summ_std.mat')
part_vals_single = part_vals_smc;

load('results_summ_continue_std.mat');
part_vals_single_continue = part_vals_smc;

nexttile

    [f,xi] = ksdensity(part_vals_summ(:,p));
    plot(xi,f,'k','LineWidth',2);
    hold on;
    [f,xi] = ksdensity(part_vals_summ_pilot(:,p));
    plot(xi,f,'b--','LineWidth',2);
    [f,xi] = ksdensity(part_vals_single(:,p));
    plot(xi,f,'g-.','LineWidth',2);
    [f,xi] = ksdensity(part_vals_single_continue(:,p));
    plot(xi,f,'r.','LineWidth',2);
plot(theta_true(p),0,'kx');

lgd = legend({'all summs','pilot (all summs)','S', 'pilot + S'});
lgd.Layout.Tile = 23;
