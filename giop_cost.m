function cost = giop_cost(guess,rin,coef)

% construct Rrs

atot = coef.aw + coef.aphs .* guess(3) + coef.adgs .* guess(1);
bbtot = coef.bbw + coef.bbps .* guess(2);
u = bbtot ./ (atot + bbtot); 
rmod = coef.g0 .* u + coef.g1 .* u .* u;

% cost function
%
cost = sum((rin-rmod).^2);
