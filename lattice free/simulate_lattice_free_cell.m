function y = simulate_lattice_free_cell(theta, extra_args)

    data = extra_args.init_data;

    if nargin < 2
        error('Not enought input arguments.');
    end
    
    if length(theta) ~= 3
        error('Incorrect number of parameter.')
    end
    if any(theta <= 0)
        error('Parameter must be positive.')
    end
    if ~isa(data, 'double') || size(data, 2) ~= 2
        error('data must be a double n by 2 matrix')
    end
    
    X = data(:, 1);
    Y = data(:, 2);
    N = size(data, 1);
    y = lattice_free_cell_mex(theta, X, Y, N, randi(1e6));
end