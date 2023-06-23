# Likelihood-Free Inference for Cosmological Parameter Inference from Weak Lensing Mass Maps 

This repository hosts the necessary code for a comprehensive project that encompasses the generation of mock weak lensing catalogues, mimicking the observed properties of the Dark Energy Survey Year 3 (DES Y3) dataset, and the measurement of various statistics from these catalogues.

main codes:

1) Mock Catalogue Generation: The repository includes code to generate mock weak lensing catalogues that closely resemble the observed properties of the DES Y3 dataset. These mock catalogues are generated based on the output of gravity-only N-body simulations, providing a realistic representation of the weak lensing mass maps for analysis (code:  populate_mocks.py)

2) Measurement of Statistics: The repository hosts code to measure different statistics from the mock catalogues. These statistics include moments, wavelet phase harmonics, scattering transform, and power spectra. By analyzing these statistics, valuable insights into the underlying cosmological parameters can be extracted. Then there are 4 main measurement codes, that measure 4 different statistics on the mocks generated with populate_mocks.py:
- run_PWH.py (wavelet phase harmonics measurement)
- run_ST.py (scattering transform measurement)
- run_cl.py (power spectrum measurement)
- run_old_moments.py (second & third moments measurement)

3) Statistics compression & posterior estimates: the code use neural compression to compress the measured statistics onto a smaller parameter space (i.e, the cosmological parameters we want to measure). It proceeds estimating the parameter posteriors of a target measurement using NDE (i.e., likelihood free inference approach). An estimate of the p-value in the compressed space is also provided. The code is in the demo_compression_LFI.ipynb notebook

# Contributing
Contributions to this project are welcome! If you encounter any issues, have ideas for improvements, or would like to add new features, please submit a pull request or open an issue on the GitHub repository.
