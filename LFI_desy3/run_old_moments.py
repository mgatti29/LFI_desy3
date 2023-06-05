import glob
import os
from Moments_analysis import moments_map
import numpy as np
import gc
import pickle
import healpy as hp


def save_obj(name, obj):
    """
    Saves an object to a pickle file.

    Args:
        name (str): The name of the pickle file.
        obj (object): The object to be saved.
    """
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, protocol=2)
        f.close()

def load_obj(name):
    """
    Loads an object from a pickle file.

    Args:
        name (str): The name of the pickle file.

    Returns:
        object: The loaded object.
    """
    with open(name + '.pkl', 'rb') as f:
        mute = pickle.load(f)
        f.close()
    return mute

def g2k_sphere(gamma1, gamma2, mask, nside=1024, lmax=2048, nosh=True):
    """
    Converts shear to convergence on a sphere using HEALPix maps.

    Args:
        gamma1 (array): Gamma1 shear map.
        gamma2 (array): Gamma2 shear map.
        mask (array): Mask for the maps.
        nside (int, optional): HEALPix nside parameter. Defaults to 1024.
        lmax (int, optional): Maximum multipole moment. Defaults to 2048.
        nosh (bool, optional): Flag indicating whether to remove the monopole and dipole terms. Defaults to True.

    Returns:
        tuple: E-mode convergence map, B-mode convergence map, and E-mode alms.
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



def compute_phmoments(file,output=''):
      
        #label = 
        #dict_temp = load_obj(file.split('pkl')[0])
        
        dict_temp = np.load(file,allow_pickle=True).item()
# 

        target = file.split('/pscratch/sd/m/mgatti/Dirac/')[1].split('.npy')[0]

            
            
        mock = target.split('runs')[1].split('_')[0]
        run = int(target.split('run')[2].split('_')[0])
        f = open(('/global/u2/m/mgatti/Mass_Mapping/peaks/params_run_1_Niall_{0}.txt'.format(mock)),'r')
        om_ = []
        h_ = []
        ns_ = []
        ob_ = []
        s8_ = []
        w_ = []
        for i,f_ in enumerate(f):
            if i>0:
                om_.append(float(f_.split(',')[0]))
                h_.append(float(f_.split(',')[4]))
                ns_.append(float(f_.split(',')[5]))
                ob_.append(float(f_.split(',')[3]))
                s8_.append(float(f_.split(',')[1]))
                w_.append(float(f_.split(',')[2]))
            else:
                print (f_)


        for rel in range(4):
            params = dict()
            params['om'] = om_[run-1]
            params['h'] = h_[run-1]
            params['s8'] = s8_[run-1]
            params['w'] = w_[run-1]
            params['ob'] = ob_[run-1]
            params['ns'] = ns_[run-1]
            for k in dict_temp[rel]['nuisance_params'].keys():
                params[k] = dict_temp[rel]['nuisance_params'][k]

            
            if not os.path.exists(output+target+'_rel{0}'.format(rel)+'.pkl'):

                conf = dict()
                conf['J_min'] = 5
                conf['J'] = 6
                conf['B'] = 2
                conf['L'] = 2
                conf['nside'] = 512
                conf['lmax'] = conf['nside']*2
                conf['verbose'] = False
                conf['output_folder'] = output_intermediate+'/test_'+target+'_rel{0}'.format(rel)

                if not os.path.exists(conf['output_folder']):
                    os.mkdir(conf['output_folder'])
                if not os.path.exists(conf['output_folder']+'/smoothed_maps'):
                    os.mkdir(conf['output_folder']+'/smoothed_maps')


    
                mcal_moments = moments_map(conf)
                #mask_DES_y3 = load_obj('/global/cfs/cdirs/des//mass_maps/Maps_final//mask_DES_y3')
                mcal_moments.conf['smoothing_scales'] = np.array([8.2,13.1,21.0,33.6,54.,86.,138,221.])
                #mask_DES_y3
                # this add the maps
                tomo_bins = [0,1,2,3] #0,1,2,3]
                for t in tomo_bins:
                    
                    e1 = np.zeros(hp.nside2npix(conf['nside']))
                    e2 = np.zeros(hp.nside2npix(conf['nside']))
                    e1n = np.zeros(hp.nside2npix(conf['nside']))
                    e2n = np.zeros(hp.nside2npix(conf['nside']))
                    e1[dict_temp[rel][t+1]['pix']] = dict_temp[rel][t+1]['e1']
                    e2[dict_temp[rel][t+1]['pix']] = dict_temp[rel][t+1]['e2']
                    e1n[dict_temp[rel][t+1]['pix']] = dict_temp[rel][t+1]['e1n']
                    e2n[dict_temp[rel][t+1]['pix']] = dict_temp[rel][t+1]['e2n']
                    
                    mask_sims = np.in1d(np.arange(len(e1)),dict_temp[rel][t+1]['pix'])
                    if rel>1:
                        f,fb,almsE    =  g2k_sphere(-e1,-e2, mask_sims, nside=conf['nside'], lmax=conf['nside']*2 ,nosh=True)
                        fn,fbn, almsEN   =  g2k_sphere(-e1n,-e2n, mask_sims, nside=conf['nside'], lmax=conf['nside']*2 ,nosh=True)
                    else:
                        f,fb,almsE    =  g2k_sphere(e1,e2, mask_sims, nside=conf['nside'], lmax=conf['nside']*2 ,nosh=True)
                        fn,fbn, almsEN   =  g2k_sphere(e1n,e2n, mask_sims, nside=conf['nside'], lmax=conf['nside']*2 ,nosh=True)

                    #fn[~mask_sims] = 0
                    
                    #sources_cat[rot][tomo_bin] = {'kE':EE,'kE_noise':EEn,'mask':mask_sims}

                    
                
                    #f = dict_temp[rel][t+1]['kE']
                    #f[~dict_temp[rel][t+1]['mask']] = 0.
                    #
                    #fn = dict_temp[rel][t+1]['kE_noise']
                    #fn[~dict_temp[rel][t+1]['mask']] = 0.
                    
                    mcal_moments.add_map(f, field_label = 'k', tomo_bin = t)
                 
                    mcal_moments.add_map(fn, field_label = 'kn', tomo_bin = t)
                    mcal_moments.add_map(fb, field_label = 'bk', tomo_bin = t)
                    mcal_moments.add_map(fbn, field_label = 'bkn', tomo_bin = t)
                    
                    
                    if t == 3:
                        mcal_moments.mask = mask_sims 
                       # mcal_moments.mask = dict_temp[rel][t+1]['mask'] 

    
                mcal_moments.transform_and_smooth('k_sm','k', shear = False, tomo_bins = [0,1,2,3], overwrite = False, skip_loading_smoothed_maps = True) 
                mcal_moments.transform_and_smooth('bk_sm','bk', shear = False, tomo_bins = [0,1,2,3], overwrite = False, skip_loading_smoothed_maps = True)   
                mcal_moments.transform_and_smooth('kn_sm','kn', shear = False, tomo_bins = [0,1,2,3], overwrite = False, skip_loading_smoothed_maps = True)   
                mcal_moments.transform_and_smooth('bkn_sm','bkn', shear = False, tomo_bins = [0,1,2,3], overwrite = False, skip_loading_smoothed_maps = True)   
            
            
           
                mcal_moments.compute_moments('KK','k_sm_kE',field_label2='k_sm_kE', tomo_bins1 = [0,1,2,3], tomo_bins2 = [0,1,2,3])
                mcal_moments.compute_moments('NN','kn_sm_kE',field_label2='kn_sm_kE', tomo_bins1 = [0,1,2,3], tomo_bins2 = [0,1,2,3])
                mcal_moments.compute_moments('bKK', 'bk_sm_kE',field_label2='bk_sm_kE', tomo_bins1 = [0,1,2,3], tomo_bins2 = [0,1,2,3])
                mcal_moments.compute_moments('bNN', 'bkn_sm_kE',field_label2='bkn_sm_kE', tomo_bins1 = [0,1,2,3], tomo_bins2 = [0,1,2,3])
                
                

                del mcal_moments.smoothed_maps
                del mcal_moments.fields
                gc.collect()


                import shutil
                shutil.rmtree(output_intermediate+'/test_'+target+'_rel{0}'.format(rel))
                save_obj(output+target+'_rel{0}'.format(rel),[mcal_moments,params])
 

#srun --nodes=4 --tasks-per-node=64   python run_old_moments_2.py 


if __name__ == '__main__':

    output_intermediate = '/pscratch/sd/m/mgatti/PWHM/temp/'
    output = '//global/cfs/cdirs/des/mgatti/Dirac/output_OLD_2_moments_new_dirac_B/' 
    files = glob.glob('/pscratch/sd/m/mgatti/Dirac/*')
    f_= []
    for f in files:
        if '512' in f:
            f_.append(f)

    runstodo = []
    count =0
    for f in f_:
        target = f.split('/pscratch/sd/m/mgatti/Dirac/')[1].split('.npy')[0]
        if not os.path.exists(output+target+'_rel3.pkl'):
            runstodo.append(f)
        else:
            count +=1



    
    print (len(runstodo),count)
    run_count=0
    #compute_phmoments(runstodo[run_count],output)
    from mpi4py import MPI 
    while run_count<len(runstodo):
        comm = MPI.COMM_WORLD
#
        if (run_count+comm.rank)<len(runstodo):
            compute_phmoments(runstodo[run_count+comm.rank],output)
                
        #if (run_count)<len(runstodo):
        #    make_maps(runstodo[run_count])
        run_count+=comm.size
        comm.bcast(run_count,root = 0)
        comm.Barrier() 
#
