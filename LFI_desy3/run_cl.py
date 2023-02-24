import numpy as np
import glob
import pymaster as nmt
import harmony as harm
import os
import healpy as hp

def g2k_sphere(gamma1, gamma2, mask, nside=1024, lmax=2048,nosh=True):
    """
    Convert shear to convergence on a sphere. In put are all healpix maps.
    """

    gamma1_mask = gamma1 * mask
    gamma2_mask = gamma2 * mask

    KQU_masked_maps = [gamma1_mask, gamma1_mask, gamma2_mask]
    alms = hp.map2alm(KQU_masked_maps, lmax=lmax, pol=True)  # Spin transform!


    ell, emm = hp.Alm.getlm(lmax=lmax)
    if nosh:
        almsE = alms[1] * 1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5
        almsB = alms[2] * 1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5
    else:
        almsE = alms[1] * 1.
        almsB = alms[2] * 1. 
    almsE[ell == 0] = 0.0
    almsB[ell == 0] = 0.0
    almsE[ell == 1] = 0.0
    almsB[ell == 1] = 0.0



    almssm = [alms[0], almsE, almsB]


    kappa_map_alm = hp.alm2map(almssm[0], nside=nside, lmax=lmax, pol=False)
    E_map = hp.alm2map(almssm[1], nside=nside, lmax=lmax, pol=False)
    B_map = hp.alm2map(almssm[2], nside=nside, lmax=lmax, pol=False)

    return E_map, B_map, almsE

import healpy as hp
config = dict()
config['nside'] = 512
# apodize mask
def doit(ii,ell_eff):
    cat = np.load(files[ii],allow_pickle=True).item()
    name = files[ii].split('/pscratch/sd/m/mgatti/Dirac/')[1].split('.npy')[0]

    for rel in range(4):
        path = root+name+'_rel{0}.npy'.format(rel)
        if not os.path.exists(path):
            cls_e1e2 = dict()
            '''
            for bin1 in range(1,5):
                for bin2 in range (bin1,5):
                    binx = '{0}_{1}'.format(bin1,bin2)
                    #print (np.sum(mask))
                    #print (cat[rel][bin1]['kE'])

                    e1 = np.zeros(hp.nside2npix(config['nside']))
                    e2 = np.zeros(hp.nside2npix(config['nside']))
                    e1n = np.zeros(hp.nside2npix(config['nside']))
                    e2n = np.zeros(hp.nside2npix(config['nside']))
                    e1[cat[rel][bin1]['pix']] = cat[rel][bin1]['e1']
                    e2[cat[rel][bin1]['pix']] = cat[rel][bin1]['e2']
                    e1n[cat[rel][bin1]['pix']] = cat[rel][bin1]['e1n']
                    e2n[cat[rel][bin1]['pix']] = cat[rel][bin1]['e2n']
                    mask = np.in1d(np.arange(len(e1)),cat[rel][bin1]['pix'])
                    be1 = np.zeros(hp.nside2npix(config['nside']))
                    be2 = np.zeros(hp.nside2npix(config['nside']))
                    be1n = np.zeros(hp.nside2npix(config['nside']))
                    be2n = np.zeros(hp.nside2npix(config['nside']))
                    be1[cat[rel][bin2]['pix']] = cat[rel][bin2]['e1']
                    be2[cat[rel][bin2]['pix']] = cat[rel][bin2]['e2']
                    be1n[cat[rel][bin2]['pix']] = cat[rel][bin2]['e1n']
                    be2n[cat[rel][bin2]['pix']] = cat[rel][bin2]['e2n']      
                    
                    f,_,_   =  g2k_sphere(e1,e2, mask==mask, nside=config['nside'], lmax=config['nside']*2 ,nosh=True)
                    fn,_,_ =  g2k_sphere(e1n,e2n,  mask==mask, nside=config['nside'], lmax=config['nside']*2 ,nosh=True)
                    
                    bf,_,_   =  g2k_sphere(be1,be2,  mask==mask, nside=config['nside'], lmax=config['nside']*2 ,nosh=True)
                    bfn,_,_ =  g2k_sphere(be1n,be2n,  mask==mask, nside=config['nside'], lmax=config['nside']*2 ,nosh=True)
                    
            # compute power spectrum
            #ls_ = dict()
            #or bin1 in range(1,5):
            #   for bin2 in range (bin1,5):
            #       binx = '{0}_{1}'.format(bin1,bin2)
            #       #print (np.sum(mask))
                    f_0 = nmt.NmtField( mask==mask, [f])
                    f_2 = nmt.NmtField( mask==mask, [bf])
                    fn_0 = nmt.NmtField( mask==mask,[fn])
                    fn_2 = nmt.NmtField( mask==mask,[bfn])


                    cl_00 = nmt.compute_full_master(f_0, f_2, b)
                    cl_00n = nmt.compute_full_master(fn_0, fn_2, b)
                    cls_e1e2[binx] = [cl_00,cl_00n,ell_eff]
            #       cls_[binx] = [cl_00,cl_00n,ell_eff]
            #       
                    
            # compute from e1e2
            '''
            cls_e1e2 = dict()
            for bin1 in range(1,5):
                for bin2 in range (bin1,5):
                    binx = '{0}_{1}'.format(bin1,bin2)
                    #print (np.sum(mask))
                    #print (cat[rel][bin1]['kE'])

                    e1 = np.zeros(hp.nside2npix(config['nside']))
                    e2 = np.zeros(hp.nside2npix(config['nside']))
                    e1n = np.zeros(hp.nside2npix(config['nside']))
                    e2n = np.zeros(hp.nside2npix(config['nside']))
                    e1[cat[rel][bin1]['pix']] = cat[rel][bin1]['e1']
                    e2[cat[rel][bin1]['pix']] = cat[rel][bin1]['e2']
                    e1n[cat[rel][bin1]['pix']] = cat[rel][bin1]['e1n']
                    e2n[cat[rel][bin1]['pix']] = cat[rel][bin1]['e2n']
                    mask = np.in1d(np.arange(len(e1)),cat[rel][bin1]['pix'])
                    be1 = np.zeros(hp.nside2npix(config['nside']))
                    be2 = np.zeros(hp.nside2npix(config['nside']))
                    be1n = np.zeros(hp.nside2npix(config['nside']))
                    be2n = np.zeros(hp.nside2npix(config['nside']))
                    be1[cat[rel][bin2]['pix']] = cat[rel][bin2]['e1']
                    be2[cat[rel][bin2]['pix']] = cat[rel][bin2]['e2']
                    be1n[cat[rel][bin2]['pix']] = cat[rel][bin2]['e1n']
                    be2n[cat[rel][bin2]['pix']] = cat[rel][bin2]['e2n']      
                    
                    
                    if rel>1:
                        f_0a = nmt.NmtField(mask, [-e1,e2])
                        f_0b = nmt.NmtField(mask, [-be1,be2])
                        f_2a = nmt.NmtField(mask, [-e1n,e2n])
                        f_2b = nmt.NmtField(mask, [-be1n,be2n])              
                    else:
                        f_0a = nmt.NmtField(mask, [e1,e2])
                        f_0b = nmt.NmtField(mask, [be1,be2])
                        f_2a = nmt.NmtField(mask, [e1n,e2n])
                        f_2b = nmt.NmtField(mask, [be1n,be2n])
                    cl_22 = nmt.compute_full_master(f_0a, f_0b, b)
                    cl_22n = nmt.compute_full_master(f_2a, f_2b, b)

      
                    cls_e1e2[binx] = [cl_22,cl_22n,ell_eff]
            
            np.save(root+name+'_rel{0}'.format(rel),cls_e1e2)

if __name__ == '__main__':

    
    #srun --nodes=4 --tasks-per-node=64 --cpus-per-task=1 --cpu-bind=cores  python compute_cl.py

    
    root = '/pscratch/sd/m/mgatti/pseudo_cl_new/'
    files_ = glob.glob('/pscratch/sd/m/mgatti/Dirac/*')



    nside = 512 # this is the nside of the maps
    _nside_h = 1024 # this one is only used to determine the binning, keep it fixed
    lmin = 8
    lmax = 2*_nside_h
    b_lmax = 3*nside-1
    n_ell_bins = 32

    b = harm.utils.make_nmtbin_powspaced(_nside_h, lmin, lmax, n_ell_bins, power=0.5, verbose=True, b_lmax=b_lmax, f_ell='pixwin') # <- f_ell should be 'pixwin' for real data, but not for simulations if the nside for the simulations matches the nside used here
    ell_eff = b.get_effective_ells()

    files = []
    
    import os
    count = 0
    nn = 0
    for file in files_:
        name = file.split('/pscratch/sd/m/mgatti/Dirac/')[1].split('.npy')[0]
        xx = root+name+'_rel{0}.npy'.format(3)
        if os.path.exists(xx):

            count +=1
        else:
            files.append(file)
        
        
    print (count,len(files))
    from mpi4py import MPI 
    run_count = 0
    while run_count<len(files):
        comm = MPI.COMM_WORLD
    #
        if (run_count+comm.rank)<len(files):
        #if (run_count)<len(files):
            #try:
               # doit(run_count,ell_eff)
                doit(run_count+comm.rank,ell_eff)
           # except:
           #     pass

        run_count+=comm.size
        comm.bcast(run_count,root = 0)
        comm.Barrier() 
