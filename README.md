# taskDemandsECoG

This repository is to reproduce the analysis and results in the following manuscript
<add final version here>

Code developed by Zeeshan Qadir, Kendrick Kay, and Dora Hermes <br> 
All scripts were tested on MacStudio M1 Ultra (macOS: Tahoe 26.6.2) and using Matlab 2025b

ECoG data collected by Muge Ozker, Daniel Yoshor, and Michael S. Beauchamp at the Baylor College of Medicine, USA <br>
ECoG data (BIDS-formatted) as well as the corresponding derivatives used in the scripts can be found here: 
<https://osf.io/tw2c3>

## Prerequisites
1.	Matlab 2023a or above
2.	Download the derivatives folder from the repository mentioned above
3.	Following dependencies: <br> 
a.	https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics.git <br> 
b.	https://github.com/cvnlab/knkutils.git <br> 
c.	https://github.com/spm/spm12.git <br> 
d.	https://github.com/fieldtrip/fieldtrip.git <br> 
e.	https://chronux.org/chronuxFiles/filesReleases/chronux_2_11.zip <br> 

## Steps to reproduce the figures:
1.	Setup your path in personalDataPath.m
a.	Setup locations for input, output and freesurfer directories
b.	Add path to mnl_ieegBasics and knkutils GitHub repositories
2.	Run from scripts preprocessing_1_baseline.m --> preprocessing_2_respElec.m and then in any order scripts fig3a.m, fig3c.m, fig4.m, fig5.m, figS3.m, figS7.m, figS8.m, figS9.m.
3.	The following scripts are each separate and run independently: figS1.m, figS2.m, figS6.m

