theta = [1,0.04,5];
init_data = csvread('dataset3/data0.csv',1,0);
rng(123)
y = simulate_lattice_free_cell(theta, init_data);

save('data_simulated/data.mat');