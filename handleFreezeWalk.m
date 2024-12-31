function [out,walk_index,walk_direction] = handleFreezeWalk(spectral_envelope, freeze_frames, hop_count, walk_index, walk_direction, freeze_buffer)

    
    if hop_count <= freeze_frames
        
        out = spectral_envelope;
    else
        walk_index = walk_index + walk_direction;
        if walk_index == freeze_frames
            walk_direction = -1;
        % If at start of buffer then switch directions
        elseif walk_index == 1
            walk_direction = 1;
        end
    out = freeze_buffer(:, walk_index);
    end






end 