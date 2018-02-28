targetRets = min(Er):0.01:max(Er);

for i =1:numel(targetRets)
    [w, sigmaP(i)] = markowitzPortfolio(Er, covMat, targetRets(i));
end

plot(sigmaP, targetRets);
grid on;
axis([0 0.25 0.05 0.25])