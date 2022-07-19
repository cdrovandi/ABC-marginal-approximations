function ssx = summStat_ma2_acf2(x, extra_args)
% summary statistic for the MA(2) example using second
% autocovariance

[T,n] = size(x);
ssx = sum(x(1:(end-2),:).*x(3:end,:));

ssx = ssx./T;

end