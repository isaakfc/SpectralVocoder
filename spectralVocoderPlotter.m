% Parameters
% x              - Source sound
% y              - Spectral envelope sound
% z              - Spectral envelope modulator sound
% lifteringType  - 'Rectangular' or 'Exponential' or 'Linear'
% morphingOn     - 1 or 0 depending if you want the sound to morph between y
% and z
% phaseEffect    - 1 or 0 Offers a choice of 3 phase effects
% whitening      - 1 or 0 Whitens the source sound
% spectralFreeze - 1 or 0 Takes an average of the spectral envelope
% spectralDelay  - 1 or 0 Employs a spectral delay for the spectral
% envelope

function [out] = spectralVocoderPlotter(x, y, z, Fs, lifteringType, morphingOn, phaseEffect, whitening, spectralFreeze, spectralDelay)

% Get sampling period
Ts = 1/Fs;

% Convert to mono
x = x(:,1);
y = y(:,1);
z = z(:,1);

% Parameters
s_win = 1024;     % Window size
n1 = 256;         % Step increment
order_que = 30;   % Cut quefrency for sound 1
r = 0.99;         % Normalizing ratio for sound output
lfoFreq = 0.5;    % Lfo frequency for morphing 
morphingMode = 0; % Initialise morphing mode
phase_mode = 0;   % Initialise phase mode
lfoType = 1;      % Initialise LFO type
    
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
L = min(length(x),length(y));
L = min(length(z),L);
pend = L - s_win; 

% Normalize audio data and zero pad
x = [zeros(s_win, 1); x; zeros(s_win-mod(L,n1),1)] / max(abs(x));
y = [zeros(s_win, 1); y; zeros(s_win-mod(L,n1),1)] / max(abs(y));
z = [zeros(s_win, 1); z; zeros(s_win-mod(L,n1),1)] / max(abs(z));

% Output audio buffer
DAFx_out = zeros(L,1);

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

    % Window the grains of source and envelope
    grain_sou = x(pin+1:pin+s_win).* w1;
    grain_env = y(pin+1:pin+s_win).* w1;
    grain_mod = z(pin+1:pin+s_win).* w1;

    %{

    if pin == n1*500
        plot(grain_env);
        xlabel('Samples');
        ylabel('Amplitude');
        title('Windowed Grain', 'FontSize', 20);
        xlim([0 1024]);
        
    end

    %}


    % Take Fourier transform of source and envelope
    f_sou = fft(grain_sou);
    f_env = fft(grain_env)/hs_win;
    f_mod = fft(grain_mod)/hs_win;
    
    % Hold phase information for later
    phase_mod = angle(f_mod);

    % Compute the cepstrum of the spectral envelope
    flog = log(0.00001+abs(f_env));

    %{

    if pin == n1*500
        f = (0:s_win/2-1)*(Fs/s_win);

        plot(f, flog(1:s_win/2));
        xlabel('Frequency (Hz)');
        ylabel('Log Amplitude');
        title('Log Frequency Spectrum', 'FontSize', 20);
        xlim([0 22000]);
        
    end
    
    %}
    

    cep = ifft(flog);
    flog_mod = log(0.00001+abs(f_mod));
    cep_mod = ifft(flog_mod);

    %{

    if pin == n1*500
        plot(cep);
        xlabel('Quefrency');
        ylabel('Cepstral Amplitude');
        title('Cepstrum', 'FontSize', 20);
        xlim([0 1024]);
        ylim([-1 1]);
        
    end

    %}

    if strcmp(lifteringType, 'Rectangular')
        % Liftering to reduce the cepstral order
        cep_cut = zeros(s_win,1);
        % Rectangular liftering
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

    
    %{
    if pin == n1*500
        plot(cep_cut);
        xlabel('Quefrency');
        ylabel('Cepstral Amplitude');
        title('Cepstrum After Exponential Liftering', 'FontSize', 20);
        xlim([0 1024]);
        ylim([-1 1]);
        
    end

    %}

    
    
    % Inverse FFT to get the spectral envelope of y
    flog_cut_env = 2*real(fft(cep_cut));
    % Inverse FFT to get the spectral envelope of z
    flog_cut_mod = 2*real(fft(cep_cut_mod));

    %{

    if pin == n1*500
        f = (0:s_win/2-1)*(Fs/s_win);
        figure();
        plot(f, flog_cut_env(1:s_win/2));
        xlabel('Frequency (Hz)');
        ylabel('Log Amplitude');
        title('Frequency Morphed Envelope', 'FontSize', 20);
        xlim([0 22000]);
        
    end
    
    %}

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
        flog_cut_env = mixSpectrums(flog_cut_env,flog_cut_mod,50);
    end


    if pin == n1*500
        f = (0:s_win/2-1)*(Fs/s_win);
        figure();
        plot(f, flog_cut_env(1:s_win/2));
        xlabel('Frequency (Hz)');
        ylabel('Log Amplitude');
        title('Frequency Morphed Envelope', 'FontSize', 20);
        xlim([0 22000]);
        
    end

    % Convert back to linear frequency (these values will only be used if
    % whitening is off)
    f_env_out_y = exp(flog_cut_env);
    f_env_out_z = exp(flog_cut_mod);

    if whitening
        % Compute cepstrum for the source
        flog_sou = log(0.00001+abs(f_sou));  % Log of the absolute value of the source FFT
        cep_sou = ifft(flog_sou);  % Inverse FFT of log of source (cepstrum of source)
        cep_cut_sou = zeros(s_win,1);  
        cep_cut_sou(1:order_que) = [cep_sou(1)/2; cep_sou(2:order_que)];  % Cut and weight the source cepstrum
        flog_cut_sou = 2*real(fft(cep_cut_sou));  % FFT of the cut source cepstrum 
        % Here, the spectral envelope is computed by subtracting the cut source cepstrum from the cut envelope cepstrum
        % This operation whitens the source by removing its spectral envelope, and then applies the envelope of the filter
        f_env_out_y = exp(flog_cut_env - flog_cut_sou);  % Whitening with source for y
        f_env_out_z = exp(flog_cut_mod - flog_cut_sou);  % Whitening with source for z
        % Morph between whitened envelopes
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

    % Overlap-add the output
    DAFx_out(pin+1:pin+s_win) = DAFx_out(pin+1:pin+s_win) + grain;

    % Move to the next frame
    pin = pin + n1;

    % Update hop count
    hop_count = hop_count + 1;
    
    % Shift delay vector
    s_delay_vector = circshift(s_delay_vector, 1);

end 

    % Trim and normalize the output
    DAFx_out = DAFx_out(s_win+1:length(DAFx_out)) / max(abs(DAFx_out));
    %sound(DAFx_out, Fs); 

    % Normalize for wav output and write to file
    out = r * DAFx_out/max(abs(DAFx_out)); 
    
end