function bondPrice = bondPriceCIR(r, kappa, theta, sigma, T, NT, nSimul)

paths = zeros(nSimul, NT+1);
paths(:,1) = r;
dt = T/NT;

for simul = 1:nSimul
    for i = 1:NT
        paths(simul, i+1) = paths(simul, i) + kappa*(theta - paths(simul, i))*dt + sigma*sqrt(paths(simul, i))*sqrt(dt)*randn();
    end
end

integral_r = dt*sum(paths, 2);
bondPrice = mean(exp(-integral_r));

end