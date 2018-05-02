function [optionValue] = trinomialOptionModel(S, K, r, T, sigma, callPut, N, euroAmer, BUpper, BLower)
%trinomial tree model to value barrier (knock-out) option

%S = stock price
%K = strike price
%r = risk free rate
%T = maturity (in years)
%sigma = volatility
%callPut = 1=Call
%N = number of steps
%euroAmer = 1=american
%BUpper = upper knock-out barrier
%BLower = lower knock-out barrier

%Step size
dt = T/N;

%Discount factor
df = exp(-r*dt);

%Find parameters
[u, pu, pm, pd] = getTreeParameters(r, sigma, dt);

%Create stock tree
stockTree = getTreeStructure(S, N, u);

%Create option tree
optionTree = getOptionTree(df, pu, pm, pd, stockTree, K, callPut, N, euroAmer, BUpper, BLower);

%Option value
optionValue = optionTree(N+1,1);

end

function [u, pu, pm, pd] = getTreeParameters(r, sigma, dt)
%Calculates trinomial tree parameters
dx = sigma*sqrt(3*dt);
u = exp(dx);

m = r - 0.5*(sigma^2);

pu = 0.5*((((sigma^2)*dt+(m^2)*(dt^2))/(dx^2)) + ((m*dt)/dx));
pd = 0.5*((((sigma^2)*dt+(m^2)*(dt^2))/(dx^2)) - ((m*dt)/dx));
pm = 1 - pu - pd;
end

function stockTree = getTreeStructure(S, N, u)
%Creates trinomial tree structure

stockTree = zeros(2*N+1, N+1);

for i = 0:N %time
    for j = -i:i %price movements
        stockTree(N+1-j, i+1) = S*(u^j);
    end
end    
end

function optionTree = getOptionTree(df, pu, pm, pd, stockTree, K, callPut, N, euroAmer, BUpper, BLower)
%Creates trinomial option tree

optionTree = stockTree*0;

%Option value at maturity
optionTree(:,end) = getPayoff(stockTree(:,end), K, callPut);

for i = N-1:-1:0 % time (backwards)
    for j = -i:i %price movement
        discountedVal = df*(pd*optionTree(N+2-j, i+2) + pm*optionTree(N+1-j, i+2) + pu*optionTree(N-j, i+2));
        if euroAmer == 1 %Amer
            excerciseVal = getPayoff(stockTree(N+1-j, i+1), K, callPut);
            optionTree(N+1-j, i+1) = max(discountedVal, excerciseVal);
        else %Euro
            optionTree(N+1-j, i+1) = discountedVal;
        end
        % Check if barrier was reached
        if stockTree(N+1-j, i+1) >= BUpper
            optionTree(N+1-j, i+1) = 0;
        elseif stockTree(N+1-j, i+1) <= BLower
            optionTree(N+1-j, i+1) = 0;
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