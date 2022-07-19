function y = simulate_ma2(theta,extra_args)
    % simulate MA(2) time series model.  The function allows to simulate
    % for than one independent replicate via the reps parameter
    
    T = extra_args.T;
    reps = extra_args.reps;
    
    e = randn(T+2,reps);
    y = e(3:end,:) + theta(1)*e(2:end-1,:) + theta(2)*e(1:end-2,:);
end


