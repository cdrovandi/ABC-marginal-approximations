
rng(7);

n=10;
y = normrnd(0,1,1,n);
ybar = mean(y);

alpha = 1;
beta = 1;
phi0 = 1^2;
mu0 = 0;


% marginal for mu

mu_grid = -3:0.01:3;
dy = zeros(1,length(mu_grid));
for i = 1:length(mu_grid)
    dy(i) = normpdf(mu_grid(i), 0 , sqrt(phi0))*(beta + 0.5*sum((y-mu_grid(i)).^2))^(-alpha - n/2);
end

dyb = normpdf(mu_grid, 0 , sqrt(phi0)).*(beta + 0.5*n*(ybar-mu_grid).^2).^(-alpha - 0.5);

% approximate normalising constants numerically
Cy = 1/trapz(mu_grid, dy);
Cyb = 1/trapz(mu_grid, dyb);

figure;
plot(mu_grid, Cy*dy, 'LineWidth', 2);
hold on;
plot(mu_grid, Cyb*dyb, 'LineWidth', 2, 'LineStyle','--');
xlim([-2 2]);
legend('\pi(\mu|ybar,s^2)','\pi(\mu|ybar)');
xlabel('\mu','FontSize',16);
ylabel('density','FontSize',16)


% marginal for phi

phi_grid = 0.01:0.01:5;
dy = zeros(1,length(phi_grid));
for i = 1:length(phi_grid)
    phi = phi_grid(i);
    dy(i) = phi^(-alpha-1)*exp(-beta/phi)*phi^(-n/2)*sqrt(2*pi/(n/phi + 1/phi0))*exp(0.5*(n*ybar/phi + mu0/phi0)^2/(n/phi + 1/phi0))*exp(-0.5*sum(y.^2)/phi);
end

Cy = 1/trapz(phi_grid, dy);
alpha_post = alpha + (n-1)/2;
beta_post = beta + 0.5*sum((y-ybar).^2);

figure;
plot(phi_grid, Cy*dy, 'LineWidth', 2);
hold on;
plot(phi_grid, beta_post^alpha_post/gamma(alpha_post)*phi_grid.^(-alpha_post - 1).*exp(-beta_post./phi_grid), 'LineWidth', 2, 'LineStyle','--');
xlim([0 4]);
legend('\pi(\phi|ybar,s^2)','\pi(\phi|s^2)');
xlabel('\phi','FontSize',16);
ylabel('density','FontSize',16)

