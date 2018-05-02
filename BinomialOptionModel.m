function [optionValue, delta, gamma, theta] = BinomialOptionModel(S, K, r, T, sigma, callPut, N, modelType, euroAmer, tau, div)
%binomial tree for pricing of European and American options

%S = stock price
%K = strike price
%r = risk free rate
%T = maturity (in years)
%sigma = volatility
%callPut = 1=Call
%N = number of steps
%modelType = 1=CRR, 2JR
%euroAmer = 1=american
%tau = time of dividend
%div = dividend size (in proportion to the stock price)

%Calculate the step size
dt = T/N;
%discount factor
df = exp(-r*dt);

% Find the tree parameters
[p, u, d] = getTreeParameters(r, sigma, modelType, dt);

%Create stock price tree
[stockTree] = getTreeStructure(S, N, u, d, dt, tau, div);

%Create option tree
[optionTree] = getOptionTree(df, p, stockTree, K, callPut, N, euroAmer);

%Calculate option value
optionValue = optionTree(end,1);

%Calculate option Greeks
[delta, gamma, theta] = getOptionGreeks(optionTree, stockTree, dt);

end

function [p, u, d] = getTreeParameters(r, sigma, modelType, dt)
% Function finds the params of the binomial tree
if (modelType==1) %CRR
    u = exp(sigma*sqrt(dt));
    d = 1/u;
    p = exp(r*dt) - d/(u-d);
elseif (modelType==2) %JR
    p = 0.5;
    u = exp((r-(1/2)*sigma^2)*dt + sigma*sqrt(dt));
    d = exp((r-(1/2)*sigma^2)*dt - sigma*sqrt(dt));
end
end

function [stockTree] = getTreeStructure(S, N, u, d, dt, tau, div)
%function creates the stock tree structure

stockTree = zeros(N+1, N+1);

for i = 0:N % time
    t = i*dt;
    for j = 0:i %price movement
        stockTree(N+1-j, i+1) = S*(u^j)*(d^(i-j))*((1-div)^(t>=tau));
    end
end
end

function [optionTree] = getOptionTree(df, p, stockTree, K, callPut, N, euroAmer)
%function calculates the option price tree

optionTree = stockTree*0;

%option prices oat maturity
optionTree(:,end) = getPayoff(stockTree(:,end), K, callPut);

for i = N-1:-1:0 %time (backwards)
    for j = 0:i %prices
       discountedVal = df*(p*optionTree(N-j, i+2) + (1-p)*optionTree(N+1-j, i+2));
       if euroAmer == 1
           exerciseVal = getPayoff(stockTree(N+1-j, i+1), K, callPut);
           optionTree(N+1-j, i+1) = max(discountedVal, exerciseVal);
       else
           optionTree(N+1-j, i+1) = discountedVal;
       end
    end    
end
end

function payoff = getPayoff(S, K, callPut)
if (callPut==1)
    payoff = max(S-K, 0);
else
    payoff = max(K-S, 0);
end
end

function [delta, gamma, theta] = getOptionGreeks(optionTree, stockTree, dt)
% Calculate option greeks

delta = (optionTree(end-1, 2) - optionTree(end, 2))/(stockTree(end-1, 2) - stockTree(end, 2));

delta_upper = (optionTree(end-2, 3) - optionTree(end-1, 3))/(stockTree(end-2, 3) - stockTree(end-1, 3));
delta_lower = (optionTree(end-1, 3) - optionTree(end, 3))/(stockTree(end-1, 3) - stockTree(end, 3));
gamma = (delta_upper - delta_lower)/(0.5*(optionTree(end-2, 3)-optionTree(end, 3)));

theta = (optionTree(end-1, 3) - optionTree(end, 1))/(2*dt);

end