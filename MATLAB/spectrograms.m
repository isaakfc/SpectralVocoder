[sig, Fs] = audioread('PhaseWhisperWhitened.wav');

t = (0:(length(sig)-1))/Fs;

Nspec = 256;
wspec = hamming(Nspec);
Noverlap = Nspec/2;

[SSig, fspec, tspec] = spectrogram(sig,wspec,Noverlap,Nspec,Fs);

PhaseSSig = angle(SSig); 

imagesc(tspec, fspec, PhaseSSig); 
axis xy; 
title('Whisper Phase Spectrum', 'FontSize', 18);
xlabel('Time (sec)', 'FontSize', 16);
ylabel('Frequency (Hz)', 'FontSize', 16);
colormap('jet');
