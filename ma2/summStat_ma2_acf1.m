function ssx = summStat_ma2_acf1(x, extra_args)
% summary statistic for the MA(2) example using variance and first
% autocovariance

[T,n] = size(x);
ssx = zeros(n,2);

for k = 1:2
   ssx(:,k) = sum(x(1:(end-(k-1)),:).*x(k:end,:));
end
ssx = ssx./T;

end