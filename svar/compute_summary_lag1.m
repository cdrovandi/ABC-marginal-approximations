function S = compute_summary_lag1(Y,extra_args)

% code to compute the lag 1 correlation summary in the sparse VAR example
% given a particular pair of (Y_i, Y_j)

pairs = extra_args.pairs;
n = extra_args.n;

pos = extra_args.pos;
x = ceil(pos/2);
pair = pairs(x,:);

iseven = rem(pos,2) == 0;

if (~iseven)
    S = mean(Y(pair(1),2:n).*Y(pair(2),1:(n-1)));
else
    S = mean(Y(pair(2),2:n).*Y(pair(1),1:(n-1)));
end

end