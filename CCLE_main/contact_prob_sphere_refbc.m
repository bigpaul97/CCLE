function [cp_matrix] = contact_prob_sphere_refbc(neff_matrix,DSBZ,a,kl)
% Calculate the contact probability between pairwise loci
% Using the new probability function of constrained polymer in a sphere 
% a = sphere radius
% kl = Kuhn length

%% Parameters

lmax = 200;
pmax = 50;

cp_matrix = neff_matrix * 0.;

%% Calculation

for l = 0:lmax
    for p = 1:pmax
        cp_matrix = cp_matrix + ((2*l+1) .* exp(-DSBZ(l+1,p)^2*kl^2/6/a^2.*neff_matrix));
    end
end

cp_matrix = cp_matrix + 1; % first zero of deriv of l=0 term (the constant term)
cp_matrix = cp_matrix .* 3/(4*pi*a^3);


end