function [phase_mode] = handlePhaseOptions(phase_on, phase_mode)


mode = input('Select phase mode, 1 for phase transfer, 2 for phase randomisation and 3 for zero phase: ');
if mode == 1
    phase_mode = 1;
elseif mode == 2
    phase_mode = 2;
elseif mode == 3
    phase_mode = 3;
else
    error('INVALID INPUT');
end
