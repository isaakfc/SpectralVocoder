function [morphingMode, lfoFreq, lfoType] = handleMorphingOptions(morphingMode, lfoFreq,lfoType)

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

mode = input('Select LFO type, 1 for sine, 2 for square, 3 for sawtooth: ');
if mode == 1
    lfoType = 1;
elseif mode == 2
    lfoType = 2;
elseif mode == 3
    lfoType = 3;
else
    error('INVALID INPUT');
end








end