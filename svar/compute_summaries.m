function S = compute_summaries(Y,extra_args)

% code to compute all 21 summaries in the sparse VAR example

nvars = extra_args.nvars;
pairs = extra_args.pairs;
n = extra_args.n;

S = zeros(1,nvars+1);

for i = 1:(nvars/2)
    pair = pairs(i,:);

    S((i-1)*2+1) = mean(Y(pair(1),2:n).*Y(pair(2),1:(n-1)));
    S(i*2) = mean(Y(pair(2),2:n).*Y(pair(1),1:(n-1)));

end

S(end) = std(Y(:));


end