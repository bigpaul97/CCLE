function [cross_corr] = chipseq_correlation(target_data,aligning_data,max_shift,resln)
% Calculate the correlation between two ChIP-seq data

% target_data: chipseq data template, with "correct" register
% aligning_data: chipseq data with "incorrect" register
% optimal_shift: what is the best register shift of aligning_data that
% gives the highest cross correlation
% max_shift (positive): the maximum register shift to test in each direction
% The picture: aligning_data is moving while target_data is fixed
% negative: aligning data is shifted leftward. vice versa

% Initialize cross correlation
cross_corr = zeros(2*max_shift,1);

% get dimension
total_size = size(target_data,1);

% get the signal directly
target_data = target_data(:,1); % change to indicate which column has the data
aligning_data = aligning_data(:,1);

% Normalize to average line
target_data = target_data - mean(target_data);
aligning_data = aligning_data - mean(aligning_data);

% Main loop
register_shift = -max_shift;
while register_shift <= max_shift
    if register_shift <= 0
        cross_corr(register_shift+max_shift+1) = (mean(aligning_data(1-register_shift:end) .* target_data(1:total_size+register_shift)) - mean(aligning_data(1-register_shift:end)) * mean(target_data(1:total_size+register_shift)))/...
            (sqrt(mean(aligning_data(1-register_shift:end).^2) - mean(aligning_data(1-register_shift:end))^2) * sqrt(mean(target_data(1:total_size+register_shift).^2) - mean(target_data(1:total_size+register_shift))^2));
    end
    
    if register_shift > 0
        cross_corr(register_shift+max_shift+1) = (mean(aligning_data(1:total_size-register_shift) .* target_data(1+register_shift:end)) - mean(aligning_data(1:total_size-register_shift)) * mean(target_data(1+register_shift:end)))/...
            (sqrt(mean(aligning_data(1:total_size-register_shift).^2) - mean(aligning_data(1:total_size-register_shift))^2) * sqrt(mean(target_data(1+register_shift:end).^2) - mean(target_data(1+register_shift:end))^2));
    end    
    
    register_shift = register_shift + 1;
end


%% Plot cross correlation vs shift

figure()

% The x-unit is in kb
plot(-max_shift*resln/1000:resln/1000:max_shift*resln/1000,cross_corr,'linewidth',2,'color','#D95319')

set(gca,'fontsize',20,'linewidth',1.2)
xlabel('Genomic shift (kb)')
ylabel('Cross correlation')


end