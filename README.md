# Spectral Vocoder

## Table of Contents
- [Overview](#overview)
- [Architecture](#architecture)
- [Processing Methods](#processing-methods)
- [Usage](#usage)
- [Future Work](#future-work)
- [References](#references)
- [License](#license)

## Overview

A MATLAB implementation of a frequency domain-based vocoder effect using cepstral analysis. This project extends the traditional vocoder concept with spectral envelope manipulation techniques and phase effects. Key features include:

* Cepstral analysis for spectral envelope extraction
* Three liftering options (Rectangular, Exponential, Linear)
* Spectral whitening capability
* Advanced envelope effects (freeze, delay, morphing)
* Phase manipulation options

## Architecture

The vocoder system processes two signal paths: the **Carrier** (source) and the **Modulator** (filter). Whitening can optionally be applied to modify spectral characteristics. Below are refined flowcharts for each component and their combination.

### Modulator/Filter Path
This path extracts the spectral envelope from the modulator signal.

```plaintext
Modulator Signal → Windowing → FFT → |FFT| and log → IFFT → Cepstrum → Liftering → FFT → Modulator Envelope
```

### Source/Carrier Path
This path provides the excitation signal. Whitening can be applied to remove carrier characteristics.

**Without Whitening:**
```plaintext
Carrier Signal → Windowing → FFT
```

**With Whitening:**
```plaintext
Carrier Signal → Windowing → FFT → |FFT| and log → IFFT → Cepstrum → Liftering → FFT → Carrier Envelope → Subtract carrier envelope from original carrier spectrum → Whitened Carrier
```

### Combined Flow
The final step combines the processed carrier and modulator paths to generate the output signal.

```plaintext
Modulator Envelope → Combine with Carrier FFT Magnitudes → IFFT → Final Output
```

## Processing Methods

### Liftering Options
The system offers three approaches to cepstral coefficient modification, each providing different levels of spectral smoothing:


<img src="https://github.com/user-attachments/assets/3c644aa0-612d-4699-98b4-44cdd5a470bf" alt="liftering" width="600">


1. **Rectangular**: Sharp cutoff for complete removal of high quefrency content
2. **Exponential**: Smooth decay providing natural transition
3. **Linear**: Gradual reduction allowing some high quefrency content

The choice of liftering affects how much of the modulator's excitation characteristics are preserved versus removed.

### Spectral Envelope Effects

These effects uniquely operate on the modulator's spectral envelope before carrier combination:

**Spectral Freeze**
- Captures N frames of spectral envelopes
- Implements controlled random-walk playback
- Creates static timbral effects
  
  
<img src="https://github.com/user-attachments/assets/e4c422fe-2219-4612-a743-8acfc8a7dd3e" alt="Spectral Freeze" width="600">


**Spectral Delay**
- Frequency-dependent delay matrix
- Individual bin delay control
- Maintains FFT symmetry
  

<img src="https://github.com/user-attachments/assets/26b60cd6-e0fe-474a-9d74-4d7c0bba5cd8" alt="Spectral Delay" width="600">


**Spectral Morphing**
- Frequency or magnitude-based blending
- LFO-controlled crossfading
- Multiple modulation waveforms

<img width="850" alt="envelope_modulation" src="https://github.com/user-attachments/assets/1f32c86c-7888-4041-8c82-584b0eade23b" />
<br>

### Phase Processing

Phase effects modify the carrier signal's phase components:

**Phase Transfer**: Applies modulator phase characteristics  
**Robotisation**: Zero-phase synthesis for monotone effect  
**Whisperisation**: Randomized phase for breathy qualities

## Usage

Basic function call:
```matlab
[out] = spectralVocoder2(x, y, z, Fs, lifteringType, morphingOn, phaseEffect, whitening, spectralFreeze, spectralDelay)
```

### Parameters
- `x`: Source/carrier signal
- `y`: Primary modulator signal
- `z`: Secondary modulator signal
- `Fs`: Sample rate
- `lifteringType`: 'Rectangular', 'Exponential', or 'Linear'
- `morphingOn`: Enable morphing (0/1)
- `phaseEffect`: Phase effect type (0-3)
- `whitening`: Enable whitening (0/1)
- `spectralFreeze`: Enable freeze (0/1)
- `spectralDelay`: Enable delay (0/1)

### Example
```matlab
% Load audio files
[x, Fs] = audioread('Synth.wav');  % Source
y = audioread('VocalSound.wav');    % Modulator 1
z = audioread('VocalSound2.wav');   % Modulator 2

% Process with whitening and exponential liftering
[out] = spectralVocoder2(x, y, z, Fs, 'Exponential', 0, 0, 1, 0, 0);
```

## Future Work

* Smoother frequency morphing transitions
* Enhanced multi-effect combinations
* More sophisticated envelope extraction methods
* Additional modulation sources and LFO shapes
* Improved effect parameter control

## References

* Zölzer, U. (2011) DAFX: Digital Audio Effects
* Furui, S. (2018) Digital Speech Processing
* Reiss, J.D. and McPherson, A. (2014) Audio Effects: Theory, Implementation and Application

## License

This project is provided for educational and research purposes. Code and design may be used with proper attribution.
