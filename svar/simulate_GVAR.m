function Y = simulate_GVAR(theta,extra_args)

% code to simulate data from sparse VAR model (note, it is assumed the pairs have already been randomly generated)

nvars = extra_args.nvars;
pairs = extra_args.pairs;
n = extra_args.n;

X = -0.1*eye(nvars);

for i = 1:(nvars/2)
    pair = pairs(i,:);

    X(pair(1), pair(2)) = theta((i-1)*2+1);
    X(pair(2), pair(1)) = theta(i*2);

end


sigma = theta(end);

Y0 = normrnd(0,sigma,nvars,1);

Y = zeros(nvars,n);

for t = 1:n
    if (t == 1)
        Y(:,t) = X*Y0 + normrnd(0,sigma,nvars,1);
    else
        Y(:,t) = X*Y(:,t-1) + normrnd(0,sigma,nvars,1);
    end
end


end