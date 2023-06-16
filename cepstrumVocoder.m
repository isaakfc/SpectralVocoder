% UX_cross_synthesis_cepstrum.m [DAFXbook, 2nd ed., chapter 8]
% This function performs a cross-synthesis with cepstrum

% Clear variables and close figures
clc;
clear all;

% Read audio files
[DAFx_sou, SR] = audioread('Synth.wav'); % Source sound
DAFx_env = audioread('VocalSound.wav'); % Sound for spectral envelope
DAFx_mod = audioread('VocalSound2.wav'); % Sound for spectral envelope modulator
Ts = 1/SR;

% Convert to mono
DAFx_sou = DAFx_sou(:,1);
DAFx_env = DAFx_env(:,1);
DAFx_mod = DAFx_mod(:,1);

% Parameters
s_win = 1024; % Window size
n1 = 256; % Step increment
order_sou = 30; % Cut quefrency for sound 1
order_env = 30; % Cut quefrency for sound 2
order_mod = 30; % Cut quefrency for modulator
r = 0.99; % Normalizing ratio for sound output
lfoFreq = 0.5;

% Initializations
w1 = hanning(s_win, 'periodic'); % Analysis window
w2 = w1; % Synthesis window
hs_win = s_win/2; % Half window size

% Buffers for audio data
grain_sou = zeros(s_win,1); % Grain for source extraction
grain_env = zeros(s_win,1); % Grain for envelope extraction
grain_mod = zeros(s_win,1); % Grain for envelope modulator extraction

% Start and end index
pin = 0; 
L = min(length(DAFx_sou),length(DAFx_env));
L = min(length(DAFx_mod),L);
pend = L - s_win; 

% Normalize audio data and zero pad
DAFx_sou = [zeros(s_win, 1); DAFx_sou; zeros(s_win-mod(L,n1),1)] / max(abs(DAFx_sou));
DAFx_env = [zeros(s_win, 1); DAFx_env; zeros(s_win-mod(L,n1),1)] / max(abs(DAFx_env));
DAFx_mod = [zeros(s_win, 1); DAFx_mod; zeros(s_win-mod(L,n1),1)] / max(abs(DAFx_mod));

% Output audio buffer
DAFx_out = zeros(L,1);

% Cross synthesis
while pin<pend
    % Window the grains of source and envelope
    grain_sou = DAFx_sou(pin+1:pin+s_win).* w1;
    grain_env = DAFx_env(pin+1:pin+s_win).* w1;
    grain_mod = DAFx_mod(pin+1:pin+s_win).* w1;

    lfoVal = sin(2*pi*Ts*pin*lfoFreq);

    % Take Fourier transform of source and envelope
    f_sou = fft(grain_sou);
    f_env = fft(grain_env)/hs_win;
    f_mod = fft(grain_mod)/hs_win;

    % Extract phase from modulator
    %phaseMod = angle(f_mod);
    % Combine magnitude of drum FFT with phase of piano FFT
    %f_env = abs(grain_env).*exp(j*phaseMod);

    % Compute the cepstrum of the spectral envelope
    flog = log(0.00001+abs(f_env));
    cep = ifft(flog);

    % Liftering to reduce the cepstral order
    cep_cut = zeros(s_win,1);
    % Rectangular liftering
    cep_cut(1:order_sou) = [cep(1)/2; cep(2:order_sou)];
   
    % Linear liftering
    %cep_lifter = (0:s_win-1)/s_win;

    % Exponential decay lifter% Exponential decay lifter
    % cep_lifter = exp(-(0:s_win-1)/s_win); 

    % Apply cut 
    %cep_cut = cep.*cep_lifter';

    % Inverse FFT to get the spectral envelope
    flog_cut = 2*real(fft(cep_cut));
    f_env_out = exp(flog_cut);

    % Now same for modulator

    % Compute the cepstrum of the spectral envelope
    %flog = log(0.00001+abs(f_mod));
    %cep = ifft(flog);

    % Liftering to reduce the cepstral order
    %cep_cut = zeros(s_win,1);
    %cep_cut(1:order_mod) = [cep(1)/2; cep(2:order_mod)];

    % Inverse FFT to get the spectral envelope
    %flog_cut = 2*real(fft(cep_cut));
    %f_mod_out = exp(flog_cut);

    %f_combined_Env = f_env_out * lfoVal + f_mod_out * (1-lfoVal);

    % Combine the source with the new spectral envelope
    grain = (real(ifft(f_sou.*f_env_out))).*w2;

    % Overlap-add the output
    DAFx_out(pin+1:pin+s_win) = DAFx_out(pin+1:pin+s_win) + grain;

    % Move to the next frame
    pin = pin + n1;
end 

% Trim and normalize the output
DAFx_out = DAFx_out(s_win+1:length(DAFx_out)) / max(abs(DAFx_out));
soundsc(DAFx_out, SR); 

% Normalize for wav output and write to file
DAFx_out_norm = r * DAFx_out/max(abs(DAFx_out)); 
audiowrite('CrossCepstrum.wav', DAFx_out_norm, SR);
