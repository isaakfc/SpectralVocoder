function [output, prev_envelope] = handle_freeze(hop_count, spectralFreeze, spectral_envelope, prev_envelope, freeze_frames)

if spectralFreeze
    if mod(hop_count, freeze_frames) == 0
        prev_envelope = spectral_envelope;
    end
    output = prev_envelope;
else
    output = spectral_envelope;
end
end

