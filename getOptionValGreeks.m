function [value, greeks] = getOptionValGreeks(S, K, r, T, sigma, w, N, amer, modelName, dSigma, dR)
% w: w=1 call opce, w=-1 put opce
% amer: amer=1 americka opce, amer=0 evropska opce
% dSigma: relativní velikost zmìny sigma pro výpoèet Vega (default = 0.01)
% dR: relativní velikost zmìny r pro výpoèet Rho (default = 0.01)
% K: realizacni cena opce
if nargin < 11
    dR = 0.01;
end
if nargin < 10
    dSigma = 0.01;
end

dt = T/N;

[value, stockTree, optionTree, u, d] = getOptionVal(S, K, r, T, sigma, w, N, amer, modelName);

if nargout > 1
    delta = getDelta(stockTree, optionTree);
    gamma = getGamma(stockTree, optionTree);
    theta = getTheta(stockTree, optionTree, u, d, r, dt, delta, gamma, sigma);
    vega = getVega(S, K, r, T, sigma, w, N, amer, modelName, dSigma);
    rho = getRho(S, K, r, T, sigma, w, N, amer, modelName, dR);
    
    greeks = [delta, gamma, theta, vega, rho];

end
end

function [value, stockTree, optionTree, u, d] = getOptionVal(S, K, r, T, sigma, w, N, amer, modelName)
dt = T/N;

[p, u, d] = getParameters(modelName, r, sigma, dt);

stockTree = getStockTree(S, u, d, N);
optionTree = getOptionTree(stockTree, w, K, p, exp(-r*dt), N, amer);

value = optionTree(end, 1);
end

function optionTree = getOptionTree(stockTree, w, K, p, disc, N, amer)

optionTree = 0*stockTree;

optionTree(:, end) = getPayoff(stockTree(:,end) ,K, w);

for i = N-1:-1:0
    for j = 0:i
        optionTree(N+1- j,i +1) = max(amer*getPayoff(stockTree(N+1- j,i +1), K, w), disc*(p*optionTree(N+1-j-1,i+2) + (1-p)*optionTree(N+1-j,i+2)));
    end
end
optionTree
end

function payoff = getPayoff(ST, K, w)

payoff = max((ST-K)*w,0);

end

function stockTree = getStockTree(S, u, d, N)
% funkce vytvari strom (matici) hodnot podkladu Si,j

stockTree = zeros(N+1, N+1);

for i = 0:N % iterujeme pres vsechny casy
    for j = 0:i % iterujeme pres posuny nahoru
        stockTree(N+1- j,i +1) = S*u^j*d^(i-j); % Si,j
    end
end
stockTree
end


function [p, u, d] = getParameters(modelName, r, sigma, dt)
% funkce vraci parametry:
% p: pravdepodobnost narustu (pst. poklesu je 1-p )
% u: koeficient narustu
% d: koeficient poklesu

if strcmpi(modelName, 'CRR') % jedna se o model Cox-Ross-Rubinstein 1979
    u = exp(sigma*sqrt(dt));
    d = 1/u; % d = exp(-sigma*sqrt(dt))
    p = (exp(r*dt)-d)/(u-d);
elseif strcmpi(modelName, 'JR') % jedna se o model Jarrow-Rudd
    u = exp((r-0.5*sigma^2)*dt + sigma*sqrt(dt));
    d = exp((r-0.5*sigma^2)*dt - sigma*sqrt(dt));
    p = 0.5;
else % model nebyl identifikovan
    warning('neznamy model, jedu CRR');
    [p, u, d] = getParameters('CRR', r, sigma, dt); % rekurzivni fce
end

end

function delta = getDelta(stockTree, optionTree)
delta = (optionTree(end-1, 2)-optionTree(end, 2))/(stockTree(end-1, 2)-stockTree(end, 2));
end

function gamma = getGamma(stockTree, optionTree)
delta1 = (optionTree(end-2, 3)-optionTree(end-1, 3))/(stockTree(end-2, 3)-stockTree(end-1, 3));
delta2 = (optionTree(end-1, 3)-optionTree(end, 3))/(stockTree(end-1, 3)-stockTree(end, 3));
denom = 0.5*(stockTree(end-2, 3)-stockTree(end, 3));
gamma = (delta1-delta2)/denom;
end

function theta = getTheta(stockTree, optionTree, u, d, r, dt, delta, gamma, sigma)
if u*d == 1
    theta = (optionTree(end-1, 3)-optionTree(end, 1))/(2*dt);
else
    theta = r*optionTree(end, 1)-r*stockTree(end, 1)*delta-0.5*sigma^2*stockTree(end, 1)^2*gamma;
end
end
   
function vega = getVega(S, K, r, T, sigma, w, N, amer, modelName, dSigma)
sigmaU = sigma*(1+dSigma);
sigmaD = sigma*(1-dSigma);

valueU = getOptionVal(S, K, r, T, sigmaU, w, N, amer, modelName);
valueD = getOptionVal(S, K, r, T, sigmaD, w, N, amer, modelName);

vega = (valueU - valueD)/(2*sigma*dSigma);
end

function rho = getRho(S, K, r, T, sigma, w, N, amer, modelName, dR)
rU = r*(1+dR);
rD = r*(1-dR);
valueU = getOptionVal(S, K, rU, T, sigma, w, N, amer, modelName);
valueD = getOptionVal(S, K, rD, T, sigma, w, N, amer, modelName);

rho = (valueU - valueD)/(2*r*dR);
end