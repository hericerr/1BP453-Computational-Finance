function optionPrice = BlackScholesSimul(S, K, T, sigma, r, NT, nSimul, w)

paths = simulateGBM(S, T, sigma, r, NT, nSimul);

ST = paths(:,end);

payoff = max(w*(ST-K),0);

optionPrice = exp(-r*T)*mean(payoff);

end