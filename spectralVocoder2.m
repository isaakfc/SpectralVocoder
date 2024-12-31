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

function [out] = spectralVocoder2(x, y, z, Fs, lifteringType, morphingOn, phaseEffect, whitening, spectralFreeze, spectralDelay)

% Get sampling period
Ts = 1/Fs;

% Convert to mono
x = x(:,1);
y = y(:,1);
z = z(:,1);

% Parameters
s_win = 1024;     % Window size
n1 = 256;         % Hop size
order_que = 30;   % Quefrency cut
r = 0.99;         % Ratio for normalising
lfoFreq = 0.5;    % Lfo frequency for morphing 
morphingMode = 0; % Initialise morphing mode
phase_mode = 0;   % Initialise phase mode
lfoType = 1;      % Initialise LFO type
    
% Set conditions for FFT loop
w1 = hanning(s_win, 'periodic'); % Hanning analysis window
w2 = w1; % Synthesis window that is the same as analysis
hs_win = s_win/2; 

% Initialising buffers
grain_sou = zeros(s_win,1); 
grain_env = zeros(s_win,1); 
grain_mod = zeros(s_win,1); 

% Initialise for indexing through loop
pin = 0; 
L = min(length(x),length(y));
L = min(length(z),L);
pend = L - s_win; 

% Zero pad and normalise
x = [zeros(s_win, 1); x; zeros(s_win-mod(L,n1),1)] / max(abs(x));
y = [zeros(s_win, 1); y; zeros(s_win-mod(L,n1),1)] / max(abs(y));
z = [zeros(s_win, 1); z; zeros(s_win-mod(L,n1),1)] / max(abs(z));

% Output audio buffer
vox_out = zeros(L,1);

% Set spectral freeze size
freeze_frames = 25;

% Set variables for spectral freeze
walk_index = 0;
walk_direction = 1;
walk_index_mod = 0;
walk_direction_mod = 1;

% Index for counting number of frames 
hop_count = 1;

% Create delay vector for spectral delay
max_delay = 200;
s_delay_vector = round(linspace(0, max_delay, s_win/2 - 1));

% Create buffers for spectral delay
delay_buffer = zeros(s_win,max_delay);
delay_buffer_mod = zeros(s_win,max_delay);


% Set morphing mode if selected
if morphingOn
     [morphingMode, lfoFreq, lfoType] = handleMorphingOptions(morphingMode, lfoFreq,lfoType);
end

% Set number of freeze frames if on
if spectralFreeze
    freeze_frames = handleFreezeOptions(freeze_frames);
end

% Set phase mode
if phaseEffect
    phase_mode = handlePhaseOptions(phaseEffect,phase_mode);
end

% For freeze effect 
freeze_buffer = zeros(1024,freeze_frames);
freeze_buffer_mod = zeros(1024,freeze_frames);

% Cross synthesis
while pin<pend
    

    if lfoType == 1
        lfoVal = (1 + sin(2*pi*Ts*pin*lfoFreq)) / 2;
    end

    if lfoType == 2
        lfoVal = (1 + square(2*pi*Ts*pin*lfoFreq)) / 2;
    end

    if lfoType == 3
        lfoVal = (1 + sawtooth(2*pi*Ts*pin*lfoFreq)) / 2;
    end



    % Window the grains for all inputs
    grain_sou = x(pin+1:pin+s_win).* w1;
    grain_env = y(pin+1:pin+s_win).* w1;
    grain_mod = z(pin+1:pin+s_win).* w1;

    % Take FFT's
    f_sou = fft(grain_sou);
    f_env = fft(grain_env)/hs_win;
    f_mod = fft(grain_mod)/hs_win;
    
    % Hold phase information for later
    phase_mod = angle(f_mod);

    % Compute cepstrum
    flog = log(0.00001+abs(f_env));
    cep = ifft(flog);
    flog_mod = log(0.00001+abs(f_mod));
    cep_mod = ifft(flog_mod);

    if strcmp(lifteringType, 'Rectangular')
        % Rectangular liftering
        cep_cut = zeros(s_win,1);
        cep_cut(1:order_que) = [cep(1)/2; cep(2:order_que)];
        % Same for modulator
        cep_cut_mod = zeros(s_win,1);
        cep_cut_mod(1:order_que) = [cep_mod(1)/2; cep_mod(2:order_que)];

    elseif strcmp(lifteringType, 'Exponential')
         % Exponential decay lifter
         cep_lifter = exp(-(0:s_win-1)/(s_win/30)); 
         % Apply cut 
         cep_cut = cep.*cep_lifter';
         cep_cut(1:1) = cep(1)/2;
         cep_cut_mod = cep_mod.*cep_lifter';
         cep_cut_mod(1:1) = cep_mod(1)/2;

    else
          %Linear liftering
          cep_lifter = 1 - (0:s_win-1)/s_win;
          cep_cut = cep.*cep_lifter';
          cep_cut(1:1) = cep(1)/2;
          cep_cut_mod = cep_mod.*cep_lifter';
          cep_cut_mod(1:1) = cep_mod(1)/2;
    end
    
    % Inverse FFT to get the spectral envelope of y
    flog_cut_env = 2*real(fft(cep_cut));
    % Inverse FFT to get the spectral envelope of z
    flog_cut_mod = 2*real(fft(cep_cut_mod));

    if spectralDelay
        delay_buffer = circularBufferWrite(flog_cut_env,delay_buffer,hop_count);
        flog_cut_env = processSpectralDelay2(flog_cut_env,delay_buffer,hop_count,s_delay_vector);
        delay_buffer_mod = circularBufferWrite(flog_cut_mod,delay_buffer_mod,hop_count);
        flog_cut_mod = processSpectralDelay2(flog_cut_mod,delay_buffer_mod,hop_count,s_delay_vector);
    end

    % Add spectral envelopes to buffer if freeze is on and not past number
    % of freeze frames
    if spectralFreeze
        freeze_buffer = handleFreezeStorage(flog_cut_env,freeze_frames, hop_count, freeze_buffer);
        [flog_cut_env,walk_index,walk_direction] = handleFreezeWalk(flog_cut_env,freeze_frames,hop_count,walk_index,walk_direction,freeze_buffer);
        freeze_buffer_mod = handleFreezeStorage(flog_cut_mod,freeze_frames, hop_count, freeze_buffer_mod);
        [flog_cut_mod,walk_index_mod,walk_direction_mod] = handleFreezeWalk(flog_cut_env,freeze_frames,hop_count,walk_index_mod,walk_direction_mod,freeze_buffer_mod);
    end

    if morphingMode == 1
        flog_cut_env = mixSpectrums(flog_cut_env,flog_cut_mod,lfoVal*100);
    end

    % Convert back to linear frequency (these values will only be used if
    % whitening is off)
    f_env_out_y = exp(flog_cut_env);
    f_env_out_z = exp(flog_cut_mod);

    if whitening
        % Compute cepstrum for the source
        flog_sou = log(0.00001+abs(f_sou));  
        cep_sou = ifft(flog_sou);  
        cep_cut_sou = zeros(s_win,1);  
        cep_cut_sou(1:order_que) = [cep_sou(1)/2; cep_sou(2:order_que)];  
        flog_cut_sou = 2*real(fft(cep_cut_sou));  
        f_env_out_y = exp(flog_cut_env - flog_cut_sou);  % Whitening with source for y
        f_env_out_z = exp(flog_cut_mod - flog_cut_sou);  % Whitening with source for z
        f_combined_env = f_env_out_y * lfoVal + f_env_out_z * (1-lfoVal);
    else
        % Morph between non whitened envelopes
        f_combined_env = f_env_out_y * lfoVal + f_env_out_z * (1-lfoVal);
    end

    if phaseEffect
        if phase_mode == 1
            f_sou = abs(f_sou).*exp(j*phase_mod);
        end
        if phase_mode == 2
             phase = [pi*rand(s_win/2+1, 1); zeros(s_win/2-1, 1)];
             phase(s_win/2+2:end) = -flip(phase(2:s_win/2));
             f_sou = abs(f_sou).*exp(j*phase);
        end
        if phase_mode == 3
            phase = zeros(s_win,1);
            f_sou = abs(f_sou).*exp(j*phase);
        end
        
    end

    if morphingMode == 2
        % ifft of morphed envelope with source
        grain = (real(ifft(f_sou.*f_combined_env))).*w2;
    else
        % Combine the source with the spectral envelope
        grain = (real(ifft(f_sou.*f_env_out_y))).*w2;
    end 

    % overlapp-add
    vox_out(pin+1:pin+s_win) = vox_out(pin+1:pin+s_win) + grain;

    % Move by hops size
    pin = pin + n1;

    % Update hop count
    hop_count = hop_count + 1;
    
    % Shift delay vector
    s_delay_vector = circshift(s_delay_vector, 1);

end 

    % Normalise the output and trim
    vox_out = vox_out(s_win+1:length(vox_out)) / max(abs(vox_out));
    sound(vox_out, Fs); 

    % Write and normalise 
    out = r * vox_out/max(abs(vox_out)); 
    
end
