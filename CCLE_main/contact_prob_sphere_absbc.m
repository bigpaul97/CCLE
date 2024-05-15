function [lookup_T] = contact_prob_sphere_absbc(array,SBZ,a,kl)
% Calculate the contact probability between pairwise loci
% Using the new probability function of constrained polymer in a sphere 
% a = sphere radius
% kl = Kuhn length

%% Parameters & initialization

lmax = 300;
pmax = 500;

num = array * 0.;
dem = array * 0.;

%% Calculation

for l = 0:lmax
    for p = 1:pmax
        num = num + ((2*l+1) .* exp(-SBZ(l+1,p)^2*kl^2/6/a^2.*array));
    end
end

for p = 1:pmax
    dem = dem + (1/p^2) .* exp(-SBZ(1,p)^2*kl^2/6/a^2.*array);
end

lookup_T = pi/(8*a^3) * num ./ dem;

end