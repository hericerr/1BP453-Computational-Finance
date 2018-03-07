function paths = simulateGBM(S, T, sigma, mu, NT, nSimul)
%nSimul poèet simulací
%NT poèet krokù

dt = T/NT;
paths = zeros(nSimul, NT+1);
paths(:,1) = S;

for simul = 1:nSimul
   for i = 1:NT
       paths(simul, i+1) = paths(simul, i) + mu*paths(simul, i)*dt+sigma*paths(simul, i)*sqrt(dt)*randn();
   end
end

end