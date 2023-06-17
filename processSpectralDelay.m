function [flog_cut_env] = processSpectralDelay(flog_cut_env,buffer,hop_count,s_delay_vector)

len = size(buffer,2);

for index = 1:length(s_delay_vector)

    if hop_count > len
        indexD = mod(hop_count-s_delay_vector(index)-1,len) + 1;
        flog_cut_env(index,1) = buffer(index,indexD);
    end


end


end