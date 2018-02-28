function [w, sigmaP] = markowitzPortfolio(Er, covMat, targetRet)

portfolioRisk = @(w)(sqrt(w(:)'*covMat*w(:)));

n = numel(Er);

A = -eye(n, n);
b = zeros(n, 1);
Aeq = [ones(1,n); Er(:)'];
beq = [1; targetRet];

[w, sigmaP] = fmincon(portfolioRisk, ones(n,1)/n, A, b, Aeq, beq);

end