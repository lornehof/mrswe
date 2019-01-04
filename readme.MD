# Magnetic Resonance Shear Wave Elastography (MRSWE)

### Overview

This repository contains code and data used to reconstruct shear wave speed elastography measurements presented in the manuscript titled **_Efficient shear wave elastography using transient acoustic radiation force excitations and magnetic resonance displacement encoding._** Manuscript was published in Magnetic Resonance in Medicine and can be found here: https://doi.org/10.1002/mrm.27647 

A brief description with animations of the MRSWE method can be found here https://lornehof.github.io/mrswe/

Algorithm was implemented in Matlab Version 2017b. Code has been tested on macOS Version 10.13.6. 

### How to Use

To generate shear wave speed maps follow these steps:
1. Make a local copy of the repository
2. Open Matlab and execute the script **mrswe.m** 
3. When prompted select the configuration file of the dataset you wish to reconstruct
4. After the calculation is completed figures will be generated and shear wave speed maps will be saved in the folder **./results/**

