function [summ] = summarise_data_robust(data)
%%
%
% Calculation of robust summary statistics
% summ(1) = q2 : the median
% summ(2) = q3 - q1 : the interquartile range
% summ(3) = (q3 + q1 - 2*q2)/(q3 - q1) : a robust measure of skewness
% summ(4) = (e7 - e5 + e3 - e1)/(q3 - q1) : a robust measure of kurtosis
%%

    p = quantile(data,[1/8,1/4,3/8,1/2,5/8,3/4,7/8]);
    e1 = p(1); q1 = p(2); e3 = p(3); q2 = p(4); 
    e5 = p(5); q3 = p(6); e7 = p(7);

    summ(1) = q2;
    summ(2) = q3 - q1;
    summ(3) = (q3 + q1 - 2*q2)/summ(2);
    summ(4) = (e7 - e5 + e3 - e1)/summ(2);
    
end