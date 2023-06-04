# LFI_desy3

This repo hosts most of the code I am using for the likelihood-free-inference project with des y3 weak lensing data.
The code folder hosts 5 main codes:
- populate_mocks.py: this code allows to generate des-y3 like weak lensing mass maps starting from PKDGRAV3 full sky lens & convergence shells.

Then there are 4 main measurement codes, that measure 4 different statistics on the mocks generated with populate_mocks.py:
- run_PWH.py (wavelet phase harmonics measurement)
- run_ST.py (scattering transform measurement)
- run_cl.py (power spectrum measurement)
- run_old_moments.py (second & third moments measurement)
