function f = prior_pdf(theta, extra_args)

f = 0;
if (sum(theta) > -1 && diff(theta) > -1 && theta(2) > -1 && theta(2) < 1)
    f = 1;
end


end