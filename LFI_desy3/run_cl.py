import numpy as np
import glob
import pymaster as nmt
import harmony as harm
import os

import healpy as hp
config = dict()
config['nside'] = 512

# Function to process the data
def doit(ii, ell_eff, weight_map):
    cat = np.load(files[ii], allow_pickle=True).item()
    name = files[ii].split('/global/cfs/cdirs/des/mgatti/Dirac_mocks/')[1].split('.npy')[0]

    for rel in range(4):
        path = root + name + '_rel{0}.npy'.format(rel)
        if not os.path.exists(path):
            cls_ = dict()

            cls_e1e2 = dict()
            for bin1 in range(1, 5):
                for bin2 in range(bin1, 5):
                    binx = '{0}_{1}'.format(bin1, bin2)

                    # Extract relevant data from the catalog
                    e1 = np.zeros(hp.nside2npix(config['nside']))
                    e2 = np.zeros(hp.nside2npix(config['nside']))
                    e1n = np.zeros(hp.nside2npix(config['nside']))
                    e2n = np.zeros(hp.nside2npix(config['nside']))
                    e1[cat[rel][bin1]['pix']] = cat[rel][bin1]['e1']
                    e2[cat[rel][bin1]['pix']] = cat[rel][bin1]['e2']
                    e1n[cat[rel][bin1]['pix']] = cat[rel][bin1]['e1n']
                    e2n[cat[rel][bin1]['pix']] = cat[rel][bin1]['e2n']
                    mask = np.in1d(np.arange(len(e1)), cat[rel][bin1]['pix'])
                    be1 = np.zeros(hp.nside2npix(config['nside']))
                    be2 = np.zeros(hp.nside2npix(config['nside']))
                    be1n = np.zeros(hp.nside2npix(config['nside']))
                    be2n = np.zeros(hp.nside2npix(config['nside']))
                    be1[cat[rel][bin2]['pix']] = cat[rel][bin2]['e1']
                    be2[cat[rel][bin2]['pix']] = cat[rel][bin2]['e2']
                    be1n[cat[rel][bin2]['pix']] = cat[rel][bin2]['e1n']
                    be2n[cat[rel][bin2]['pix']] = cat[rel][bin2]['e2n'] 

                    # Create NmtField objects
                    f_0a = nmt.NmtField(weight_map[rel][bin1], [e1, e2])
                    f_0b = nmt.NmtField(weight_map[rel][bin2], [be1, be2])
                    f_2a = nmt.NmtField(weight_map[rel][bin1], [e1n, e2n])
                    f_2b = nmt.NmtField(weight_map[rel][bin2], [be1n, be2n])

                    # Compute power spectra
                    cl_22 = nmt.compute_full_master(f_0a, f_0b, b)
                    cl_22n = nmt.compute_full_master(f_2a, f_2b, b)

                    # Store the computed power spectra
                    cls_[binx] = [cl_22, cl_22n, ell_eff]

            np.save(root + name + '_rel{0}'.format(rel), cls_)

if __name__ == '__main__':
    root = '/global/cfs/cdirs/des/mgatti/Dirac/LFI_dv/cl/'
    files_ = glob.glob('/global/cfs/cdirs/des/mgatti/Dirac_mocks/*')

    # Load weight maps
    weight_map = np.load('/global/cfs/cdirs/des/mass_maps/Maps_final/weight_maps.npy', allow_pickle=True).item()

    nside = 512  # This is the nside of the maps
    _nside_h = 512  # This one is only used to determine the binning, keep it fixed
    lmin = 8
    lmax = 3 * _nside_h
    b_lmax = 3 * nside - 1
    n_ell_bins = 28

    # Create NmtBin object for binning the power spectra
    b = harm.utils.make_nmtbin_powspaced(_nside_h, lmin, lmax, n_ell_bins, power=0.5, verbose=True, b_lmax=b_lmax,
                                         f_ell='pixwin')
    ell_eff = b.get_effective_ells()

    files = []
    count = 0
    nn = 0
    for file in files_:
        if ('_noiserel100' in file) or ('_noiserel101' in file) or ('_noiserel102' in file):
            name = file.split('/global/cfs/cdirs/des/mgatti/Dirac_mocks/')[1].split('.npy')[0]
            xx = root + name + '_rel{0}.npy'.format(3)
            if os.path.exists(xx):
                count += 1
            else:
                files.append(file)

    print(count, len(files))

    from mpi4py import MPI

    run_count = 0
    while run_count < len(files):
        comm = MPI.COMM_WORLD

        if (run_count + comm.rank) < len(files):
            try:
                doit(run_count + comm.rank, ell_eff, weight_map)
            except:
                pass

        run_count += comm.size
        comm.bcast(run_count, root=0)
        comm.Barrier()
