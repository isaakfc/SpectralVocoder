function [freeze_buffer] = handle_freeze_storage(spectral_envelope, freeze_frames, hop_count, freeze_buffer)

    if hop_count <= freeze_frames
        freeze_buffer(:, hop_count) = spectral_envelope;
    end

end