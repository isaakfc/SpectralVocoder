[x, Fs] = audioread('Synth.wav'); % Source sound
y = audioread('VocalSound.wav'); % Sound for spectral envelope
z = audioread('VocalSound2.wav'); % Sound for spectral envelope modulator

% Parameters
% x              - Source sound
% y              - Spectral envelope sound
% z              - Spectral envelope modulator sound
% lifteringType  - 'Rectangular' or 'Exponential' or 'Linear'
% morphingOn     - 1 or 0 depending if you want the sound to morph between y
% and z
% phaseTransfer  - Gives options for phase effects
% whitening      - Whitens the source sound
% spectralFreeze - 1 or 0 Freezes a selectable amount of frames from
% beginning
% spectralDelay  - 1 or 0 Employs a spectral delay for the spectral
% envelope

[out] = spectralVocoder2(x, y, z, Fs, 'Rectangular', 0, 0, 1, 0, 0);
%[out] = spectralVocoderPlotter(x, z, y, Fs, 'Rectangular', 0, 0, 1, 0, 0);

