# Deconvolution project for state space estimation

The multi-channel deconvolution project for the ECE 6397: State space estimation.

Now the project is being processed. We will update more informations in the future.

## Matlab project for reference

This project is our baseline reference. The work is done by *Amin, MD Rafiul*, and *Rose T. Faghih*. The original readme text is as bellow:

> Please see the example matlab file "`Exp_on_real_data_3_channel_4_25_2019`"

> If you are any part of the repository please cite the following paper:

> Amin, MD Rafiul, and Rose T. Faghih. "[Robust Inference of Autonomic Nervous System Activation Using Skin Conductance Measurements: A Multi-Channel Sparse System Identification Approach.][paper-robust]" IEEE Access 7 (2019): 173419-173437.

## Description for the datasets

More datasets could be downloaded here:

[PsPM-SCRV10: Skin conductance responses to loud sounds, simultanously recorded from palm, fingers and foot][link-dset]

Repository Version: 2017.02.14

This dataset includes skin conductance response (SCR) measurements recorded from 3 locations (see section 'Data files' below) for each of 26 healthy unmedicated participants (12 males and 14 females aged 24.4+/-4.9 years) in response to an auditory stimulus (single white noise bursts, 1s length; 10ms ramp; ~85dB). Participants are asked to press a foot pedal upon hearing a stimulus.

--------------------------------------------------------------------------------
This repository is curated by the PsPM team, pspm.sourceforge.net

--------------------------------------------------------------------------------
Dataset is described and used in following references:
- Dominik R. Bach & Guillaume Flandin & Karl J. Friston & Raymond J. Dolan (2010) Modelling event-related skin conductance responses. International Journal of Psychophysiology, 75, pp 349-356. DOI:10.1016/j.ijpsycho.2010.01.005

--------------------------------------------------------------------------------
File Structure:
```bash
-Readme.txt             // This document
-+Data                  // Folder Containing all data files
 |-SCRV_10_spike_XX.mat  // 22 files holding biosignal traces 
 |-SCRV_10_cogent_XX.mat // 22 files holding keyboard presses information
 |-...
```

-------------------------------------------------------------------------------
Data Files:
- Files holding biosignal traces, recorded with CED spike software.
   - event markers used for synchronization with Cogent files, marker are set at sound presentation onset (channel 1)
   - skin conductance on the thenar/hypothenar of the non-dominant hand (channel 2)
   - skin conductance on the volar middle phalanx of the dominant 2nd/3rd finger (channel 3)
   - skin conductance on the medial plantar surface of the non-dominant foot (channel 4)
   - heart beat time stamps(channel 5), never checked for integrity.
   - respiration (channel 6), never checked for integrity.
- Files holding keyboard presses, recorded with Cogent 2000 software.

The data file name is composed of: 
 - SCRV_10       The dataset name
 - spike/cogent The type of data recorded (see above)
 - XX           The subject number

Data are stored as .mat files for use with MATLAB (The MathWorks Inc., Natick, USA) in a format readable by the PsPM toolbox (pspm.sourceforge.net). All matlab files are saved in MATLAB 8.6 (R2015b) format.

## Update report

### 0.3 @ 02/14/2020

1. Add the reference project branch.

### 0.2 @ 02/08/2020

1. Finish the first edition of the draft for the introduction and literature review.
2. Add the literature review note 2.

### 0.1 @ 01/27/2020

1. Add the notes from Jin Lu.
2. Finish the abstract.
3. Fix some typos.

### 0.1-a @ 01/26/2020

1. Begin to work on the abstract.

[paper-robust]:https://ieeexplore.ieee.org/abstract/document/8917550
[link-dset]:https://zenodo.org/record/291465#.XkcDCWhKiUl