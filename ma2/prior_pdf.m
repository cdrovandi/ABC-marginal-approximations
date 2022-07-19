function f = prior_pdf(theta, extra_args)

% pdf for prior

f = 0;
if (sum(theta) > -1 && diff(theta) > -1 && theta(2) > -1 && theta(2) < 1)
    f = 1;
end


end