function [f,lx,lp] = discretize_normal_var_kurt(y,n,mu1,mu2,mu4,sigmaT)
    lwidth = y(1);
    llambda = abs(y(2)) / (1+abs(y(2)));
    lx = linspace(-lwidth*sqrt(mu2)+mu1, lwidth*sqrt(mu2)+mu1,n)';
    lp = aux.discrete_normal_alt(lx,mu1,sigmaT);
    lmass0 = zeros(n,1);
    lmass0((n+1)/2) = 1; 
    
    lp = (1-llambda).*lp + lmass0.*llambda;

    Ex2 = sum(lp.*(lx.^2));
    Ex4 = sum(lp.*(lx.^4));
    
    f = [sqrt(Ex2)-sqrt(mu2); Ex4.^0.25 - mu4.^0.25];
    
end 

