function f = prior_sampler(extra_args)

% function to sample prior

while(1)
    f = unifrnd([-2 -1], [2 1]);
    if (sum(f) > -1 && diff(f) > -1 && f(2) > -1 && f(2) < 1)
        return;
    end
end

end