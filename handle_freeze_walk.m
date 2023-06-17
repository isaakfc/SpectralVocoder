function [out,walk_index,walk_direction] = handle_freeze_walk(spectral_envelope, freeze_frames, hop_count, walk_index, walk_direction, freeze_buffer)

    
    if hop_count <= freeze_frames
        
        out = spectral_envelope;
    else
        walk_index = walk_index + walk_direction;
        if walk_index == freeze_frames
            walk_direction = -1;
        % If we've reached the start of the buffer, change direction
        elseif walk_index == 1
            walk_direction = 1;
        end
    out = freeze_buffer(:, walk_index);
    end






end 