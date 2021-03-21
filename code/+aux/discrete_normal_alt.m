function [p] = discrete_normal_alt(x,mu,sigma)
% Creates equally spaced approximation to normal distribution
% x is the grid
% mu is mean
% sigma is standard deviation
% f is the error in the approximation
% x gives the location of the points
% p is probabilities

n = length(x);
if n==2
    p = 0.5.*ones(n,1);
elseif n>2    
    p  = zeros(n,1);
    p(1) = normcdf(x(1) + 0.5*(x(2)-x(1)),mu,sigma);
    for i = 2:n-1
        p(i) = normcdf(x(i) + 0.5*(x(i+1)-x(i)),mu,sigma) - normcdf(x(i) - 0.5*(x(i)-x(i-1)),mu,sigma);
    end
    p(n) = 1 - sum(p(1:n-1));
end

end