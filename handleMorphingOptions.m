function [morphingMode, lfoFreq] = handleMorphingOptions(morphingMode, lfoFreq)

lfo = input('Enter value for morphing LFO between 0.1 and 1: ');

if lfo < 0.1 || lfo > 1
    error('INVALID INPUT');
else
    lfoFreq = lfo;
end


mode = input('Select morphing mode, 1 for frequency morphing or 2 for magnitude morphing: ');
if mode == 1
    morphingMode = 1;
elseif mode == 2
    morphingMode = 2;
else
    error('INVALID INPUT');
end









end