import glob
import os
from Moments_analysis import moments_map
import numpy as np
import gc
import pickle
import healpy as hp
def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, protocol=2)
        f.close()

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        mute =  pickle.load(f)
        f.close()
    return mute


import sys
sys.path.append('/global/u2/m/mgatti/Mass_Mapping/peaks/scattering_transform')
import scattering

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

    return E_map, B_map, almsE, almsB

def compute_phmoments(file,output=''):
      
        #label = 
        #dict_temp = load_obj(file.split('pkl')[0])
        
        dict_temp = np.load(file,allow_pickle=True).item()
# 

        target = file.split('/global/cfs/cdirs/des/mgatti/Dirac_mocks/')[1].split('.npy')[0]

            
            
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
                conf['j_min'] = 0
                conf['J'] = 6
                conf['B'] = 2
                conf['L'] = 2
                conf['nside'] = 512
                conf['lmax'] = conf['nside']*2
                conf['verbose'] = False
                conf['output_folder'] = output_intermediate+'/test_'+target+'_rel{0}'.format(rel)


                mcal_moments = moments_map(conf)
                #mask_DES_y3 = load_obj('/global/cfs/cdirs/des//mass_maps/Maps_final//mask_DES_y3')

                #mask_DES_y3
                # this add the maps
                tomo_bins = [0,1,2,3] #0,1,2,3]
                almsB_ = []
                almsE_ = []
                almsBN_ = []
                almsEN_ = []
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
                    #if rel>1:
                    if 1==1:
                    #    f,fb,almsE , almsB   =  g2k_sphere(-e1,-e2, mask_sims, nside=conf['nside'], lmax=conf['nside']*2 ,nosh=True)
                    #    fn,fbn, almsEN , almsBN  =  g2k_sphere(-e1n,-e2n, mask_sims, nside=conf['nside'], lmax=conf['nside']*2 ,nosh=True)
                    #else:
                        f,fb,almsE , almsB   =  g2k_sphere(e1,e2, mask_sims, nside=conf['nside'], lmax=conf['nside']*2 ,nosh=True)
                        fn,fbn, almsEN , almsBN  =  g2k_sphere(e1n,e2n, mask_sims, nside=conf['nside'], lmax=conf['nside']*2 ,nosh=True)

                    almsB_.append(almsB)
                    almsE_.append(almsE)
                    almsBN_.append(almsBN)
                    almsEN_.append(almsEN)
                    #f[~mask_sims] = 0
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
      
                for t1 in tomo_bins:
                    for t2 in tomo_bins:
                        if t1>t2:
                            kE = hp.alm2map(almsE_[t1]*almsE_[t2], nside=conf['nside'], lmax=conf['nside']*2, pol=False)
                            kB = hp.alm2map(almsB_[t1]*almsB_[t2], nside=conf['nside'], lmax=conf['nside']*2, pol=False)
                            kEN = hp.alm2map(almsEN_[t1]*almsEN_[t2], nside=conf['nside'], lmax=conf['nside']*2, pol=False)
                            kBN = hp.alm2map(almsBN_[t1]*almsBN_[t2], nside=conf['nside'], lmax=conf['nside']*2, pol=False)
                            mcal_moments.add_map(kE, field_label = 'k', tomo_bin = 10*t1+t2)
                            mcal_moments.add_map(kB, field_label = 'bk', tomo_bin = 10*t1+t2)
                            mcal_moments.add_map(kEN, field_label = 'kn', tomo_bin = 10*t1+t2)
                            mcal_moments.add_map(kBN, field_label = 'bkn', tomo_bin = 10*t1+t2)


                mcal_moments.cut_patches( nside=512, nside_small=8)
                
                del mcal_moments.fields
                gc.collect()
                mcal_moments.moments_ST = dict()
                
                st_calc = scattering.Scattering2d(M=128, N=128, J=6, L=3)
                mcal_moments.moments_ST['K'] = dict()
                mcal_moments.moments_ST['BK'] = dict()
                mcal_moments.moments_ST['KN'] = dict()
                mcal_moments.moments_ST['BKN'] = dict()
           
                t_ = []
                for i in tomo_bins:
                    t_.append(i)
                for t1 in tomo_bins:
                    for t2 in tomo_bins:
                        if t1>t2:
                            t_.append(10*t1+t2)
                            
                for i in t_:
          
                    print (i)
                    s_mean = st_calc.scattering_coef(np.load(mcal_moments.fields_patches['k'][i]+'.npy',allow_pickle=True))
                    label = 'K'
                    mcal_moments.moments_ST[label][i]= dict()
                    mcal_moments.moments_ST[label][i]['S0'] = s_mean['S0'].cpu().data.numpy().mean(axis=0)
                    mcal_moments.moments_ST[label][i]['S1'] = s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten() [s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten() == s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten()]
                    mcal_moments.moments_ST[label][i]['S2'] = s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten() [s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten() == s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten()]
                    
                    #s_mean = st_calc.scattering_coef(np.load(mcal_moments.fields_patches['bk'][i]+'.npy',allow_pickle=True))
                    #label = 'BK'
                    #mcal_moments.moments_ST[label][i]= dict()
                    #mcal_moments.moments_ST[label][i]['S0'] = s_mean['S0'].cpu().data.numpy().mean(axis=0)
                    #mcal_moments.moments_ST[label][i]['S1'] = s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten() [s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten() == s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten()]
                    #mcal_moments.moments_ST[label][i]['S2'] = s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten() [s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten() == s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten()]
                    #             
                    s_mean = st_calc.scattering_coef(np.load(mcal_moments.fields_patches['kn'][i]+'.npy',allow_pickle=True))
                    label = 'KN'
                    mcal_moments.moments_ST[label][i] = dict()
                    mcal_moments.moments_ST[label][i]['S0'] = s_mean['S0'].cpu().data.numpy().mean(axis=0)
                    mcal_moments.moments_ST[label][i]['S1'] = s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten() [s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten() == s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten()]
                    mcal_moments.moments_ST[label][i]['S2'] = s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten() [s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten() == s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten()]
                                   
                    #s_mean = st_calc.scattering_coef(np.load(mcal_moments.fields_patches['bkn'][i]+'.npy',allow_pickle=True))
                    #label = 'BKN'
                    #mcal_moments.moments_ST[label][i] = dict()
                    #mcal_moments.moments_ST[label][i]['S0'] = s_mean['S0'].cpu().data.numpy().mean(axis=0)
                    #mcal_moments.moments_ST[label][i]['S1'] = s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten() [s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten() == s_mean['S1_iso'].cpu().data.numpy().mean(axis=0).flatten()]
                    #mcal_moments.moments_ST[label][i]['S2'] = s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten() [s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten() == s_mean['S2_iso'].cpu().data.numpy().mean(axis=0).flatten()]
                            
                #try:
                #    #del mcal_moments.fields
                #    del mcal_moments.fields_patches
                #    #del mcal_moments.smoothed_maps
                #    gc.collect()
                #except:
                #    pass
                #print ('save')

                #import shutil
                #shutil.rmtree(output_intermediate+'/test_'+target+'_rel{0}'.format(rel))
                save_obj(output+target+'_rel{0}'.format(rel),[mcal_moments,params])
 
                #try:
                #    #del mcal_moments.fields
                #    del mcal_moments.fields_patches
                #    #del mcal_moments.smoothed_maps
                #    gc.collect()
                #except:
                #    pass
                #print ('save')
#
                #import shutil
                #shutil.rmtree(output_intermediate+'/test_'+target+'_rel{0}'.format(rel))
                #save_obj(output+'moments_'+target+'_rel{0}'.format(rel),[mcal_moments,params])
 
#srun --nodes=4 --tasks-per-node=64   python run_ST.py 
#srun --nodes=4 --tasks-per-node=16   python run_ST.py 


if __name__ == '__main__':

    output_intermediate = '/pscratch/sd/m/mgatti/PWHM/temp/'
    output = '/global/cfs/cdirs/des/mgatti/Dirac/LFI_dv/ST/' 
    files = glob.glob('/global/cfs/cdirs/des/mgatti/Dirac_mocks/*')
    f_= []
    for f in files:
        try:
            xx = np.int(f.split('noiserel')[1].split('.npy')[0])

            if ('512' in f) :
                f_.append(f)
        except:
            pass

    runstodo = []
    count = 0
    for f in f_:
            target = f.split('/global/cfs/cdirs/des/mgatti/Dirac_mocks/')[1].split('.npy')[0]
        #if ('_noiserel6' in f) or ('_noiserel7' in f) or ('_noiserel8' in f):
        
            if not os.path.exists(output+target+'_rel3.pkl'):
                # if ('noiserel31' in target) or ('noiserel32' in target):
                    runstodo.append(f)
            else:
                count +=1



            
    #compute_phmoments(runstodo[0],output)
 


   
    print (len(runstodo),count)
    run_count=0
    from mpi4py import MPI 
    while run_count<len(runstodo):
        comm = MPI.COMM_WORLD
#
        if (run_count+comm.rank)<len(runstodo):
            try:
                compute_phmoments(runstodo[run_count+comm.rank],output)
            except:
                pass
        #if (run_count)<len(runstodo):
        #    make_maps(runstodo[run_count])
        run_count+=comm.size
        comm.bcast(run_count,root = 0)
        comm.Barrier() 
 #
