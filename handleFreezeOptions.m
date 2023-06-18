function [freeze_frames] = handleFreezeOptions(freeze_frames)

freeze_frames = input('Enter number between 10 and 50 for number of freeze frames to use: ');

if freeze_frames < 10 || freeze_frames > 50
    error('INVALID INPUT');
end



end