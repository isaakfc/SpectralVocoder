function [flog_cut_env] = processSpectralDelay2(flog_cut_env,buffer,hop_count,s_delay_vector)

len = size(buffer,2);

% Get whole length of frequency bins
bin_length = length(s_delay_vector) * 2 + 2;

for index = 1:length(s_delay_vector)

    if hop_count > len
        indexD = mod(hop_count-s_delay_vector(index)-1,len) + 1;
        % For first half of bins
        flog_cut_env(index+1,1) = buffer(index,indexD);
        % For seond half of bins
        flog_cut_env(bin_length - index + 1,1) = buffer(index,indexD);
    end


end


end