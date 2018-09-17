function optionPrice = getEurobarierGigitalValue(S, T, sigma, r, K, H, p, NT, nSimul)

paths = simulateGBM(S, T, sigma, r, NT, nSimul);
ST = paths(:,end);
payoff = zeros(nSimul,1);
for i = 1:nSimul
    if and(ST(i) >= K, ST(i) <= H)
        payoff(i) = p;
    else
        payoff(i) = 0;
    end
end

optionPrice = exp(-r*T)*mean(payoff);
end

