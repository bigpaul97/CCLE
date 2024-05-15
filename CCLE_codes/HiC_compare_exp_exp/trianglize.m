function result = trianglize(experiment,simulation,L)
% make the heatmap upper triangle to be exp. result while lower triangle to
% be simulation

result = zeros(L,L);

for i = 1:L
    for j = i+1:L
        result(i,j) = experiment(j,i);
        result(j,i) = simulation(j,i); % diagonal is zero
    end
end

end