
function [output] = handle_freeze(buffer, hop_count, spectralFreeze,spectral_envelope)

if spectralFreeze

    buffer = circularBuffer(spectral_envelope,buffer,hop_count);

    if hop_count <= size(buffer, 2)
        output = mean(buffer(:, 1:hop_count), 2);
    else
        output = mean(buffer, 2);
    end

else
    output = spectral_envelope;
end

    

end