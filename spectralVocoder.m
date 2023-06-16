% Parameters
% x              - Source sound
% y              - Spectral envelope sound
% z              - Spectral envelope modulator sound
% lifteringType  - 'Rectangular' or 'Exponential' or 'Linear'
% morphingOn     - 1 or 0 depending if you want the sound to morph between y
% and z
% phaseTransfer  - Transfers the phase of the modulator into y
% whitening      - Whitens the source sound


function [out] = spectralVocoder(x, y, z, Fs, lifteringType, morphingOn, phaseTransfer, whitening)

% Get sampling period
Ts = 1/Fs;

% Convert to mono
x = x(:,1);
y = y(:,1);
z = z(:,1);

% Parameters
s_win = 1024;   % Window size
n1 = 256;       % Step increment
order_que = 30; % Cut quefrency for sound 1
r = 0.99;       % Normalizing ratio for sound output
lfoFreq = 0.5;  % Lfo frequency for morphing 
    
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

% Cross synthesis
while pin<pend
    % Window the grains of source and envelope
    grain_sou = x(pin+1:pin+s_win).* w1;
    grain_env = y(pin+1:pin+s_win).* w1;
    grain_mod = z(pin+1:pin+s_win).* w1;

    % Take Fourier transform of source and envelope
    f_sou = fft(grain_sou);
    f_env = fft(grain_env)/hs_win;
    f_mod = fft(grain_mod)/hs_win;

    if phaseTransfer
        % Extract phase from modulator
        phase_mod = angle(f_mod);
        % Combine magnitude of drum FFT with phase of piano FFT
        f_env = abs(grain_env).*exp(j*phase_mod);
    end


    % Compute the cepstrum of the spectral envelope
    flog = log(0.00001+abs(f_env));
    cep = ifft(flog);
    flog_mod = log(0.00001+abs(f_mod));
    cep_mod = ifft(flog_mod);

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
         cep_lifter = exp(-(0:s_win-1)/s_win); 
         % Apply cut 
         cep_cut = cep.*cep_lifter';
         cep_cut_mod = cep_mod.*cep_lifter';

    else
          %Linear liftering
          cep_lifter = (0:s_win-1)/s_win;
          cep_cut = cep.*cep_lifter';
          cep_cut_mod = cep_mod.*cep_lifter';
    end
    
    % Inverse FFT to get the spectral envelope
    flog_cut_env = 2*real(fft(cep_cut));
    f_env_out = exp(flog_cut_env);

    if whitening
        % Calculate the cepstrum for the source
        % This is done by taking the inverse FFT of the log of the absolute value of the spectrum
        flog_sou = log(0.00001+abs(f_sou));  % Log of the absolute value of the source FFT
        cep_sou = ifft(flog_sou);  % Inverse FFT of log of source (cepstrum of source)
        cep_cut_sou = zeros(s_win,1);  
        cep_cut_sou(1:order_que) = [cep_sou(1)/2; cep_sou(2:order_que)];  % Cut and weight the source cepstrum
        flog_cut_sou = 2*real(fft(cep_cut_sou));  % FFT of the cut source cepstrum 

        % Here, the spectral envelope is computed by subtracting the cut source cepstrum from the cut envelope cepstrum
        % This operation whitens the source by removing its spectral envelope, and then applies the envelope of the filter
        f_env_out = exp(flog_cut_env - flog_cut_sou);  % Whitening with source

    end

    if morphingOn
        % Compute lfo value at current frame
        lfoVal = sin(2*pi*Ts*pin*lfoFreq);
        % Inverse FFT to get the spectral envelope
        flog_cut_mod = 2*real(fft(cep_cut_mod));
        f_mod_out = exp(flog_cut_mod);
        f_combined_env = f_env_out * lfoVal + f_mod_out * (1-lfoVal);
        % Combine the source with the new morphed spectral envelope
        grain = (real(ifft(f_sou.*f_combined_env))).*w2;
    else
        % Combine the source with the spectral envelope
         grain = (real(ifft(f_sou.*f_env_out))).*w2;
    end 

    % Overlap-add the output
    DAFx_out(pin+1:pin+s_win) = DAFx_out(pin+1:pin+s_win) + grain;

    % Move to the next frame
    pin = pin + n1;
end 

    % Trim and normalize the output
    DAFx_out = DAFx_out(s_win+1:length(DAFx_out)) / max(abs(DAFx_out));
    sound(DAFx_out, Fs); 

    % Normalize for wav output and write to file
    out = r * DAFx_out/max(abs(DAFx_out)); 
    
end



