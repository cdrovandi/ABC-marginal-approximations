function ssx = summStat_ma2_acf(x, extra_args)
% summary statistic for the MA(2) example based on autocovariances
% NOTE that K = 1 represents the variance

K = extra_args.K;

[T,n] = size(x);
ssx = zeros(n,K);

for k = 1:K
   ssx(:,k) = sum(x(1:(end-(k-1)),:).*x(k:end,:));
end
ssx = ssx./T;

end