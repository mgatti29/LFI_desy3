from bornraytrace import lensing
import bornraytrace
import gc
import pyfits as pf
import pickle
import numpy as np
import healpy as hp
import os
import copy
import scipy
from scipy.interpolate import interp1d
import timeit
from bornraytrace import lensing as brk
from astropy import units as u
import gc
from bornraytrace import intrinsic_alignments as iaa
from Moments_analysis import convert_to_pix_coord, IndexToDeclRa, apply_random_rotation, addSourceEllipticity, gk_inv
import frogress
from astropy.table import Table  
import pandas as pd
from astropy.cosmology import z_at_value
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
sys.path.insert(0, "/global/homes/m/mgatti/PKDGRAV/lfi_project/scripts/")
from utility import *

def random_draw_ell_from_w(wi,w,e1,e2):
    '''
    wi: input weights
    w,e1,e2: all the weights and galaxy ellipticities of the catalog.
    e1_,e2_: output ellipticities drawn from w,e1,e2.
    '''


    ell_cont = dict()
    for w_ in np.unique(w):
        mask_ = w == w_
        w__ = np.int(w_*10000)
        ell_cont[w__] = [e1[mask_],e2[mask_]]

    e1_ = np.zeros(len(wi))
    e2_ = np.zeros(len(wi))


    for w_ in np.unique(wi):
        mask_ = (wi*10000).astype(np.int) == np.int(w_*10000)
        e1_[mask_] = ell_cont[np.int(w_*10000)][0][np.random.randint(0,len(ell_cont[np.int(w_*10000)][0]),len(e1_[mask_]))]
        e2_[mask_] = ell_cont[np.int(w_*10000)][1][np.random.randint(0,len(ell_cont[np.int(w_*10000)][0]),len(e1_[mask_]))]

    return e1_,e2_


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


def rotate_map_approx(mask, rot_angles, flip=False,nside = 1024):
    alpha, delta = hp.pix2ang(nside, np.arange(len(mask)))

    rot = hp.rotator.Rotator(rot=rot_angles, deg=True)
    rot_alpha, rot_delta = rot(alpha, delta)
    if not flip:
        rot_i = hp.ang2pix(nside, rot_alpha, rot_delta)
    else:
        rot_i = hp.ang2pix(nside, np.pi-rot_alpha, rot_delta)
    rot_map = mask*0.
    rot_map[rot_i] =  mask[np.arange(len(mask))]
    return rot_map



'''
This code takes the output of pkdgrav sims. It generates simulated des y3 like catalogs, adding shape noise and weights from the fiducial des y3 catalog on data. 

how to run it in parallel (when you have multiple  sims:
srun --nodes=4 --tasks-per-node=32 --cpus-per-task=2 --cpu-bind=cores  python populate_mocks.py



srun --nodes=4 --tasks-per-node=10 --cpus-per-task=6 --cpu-bind=cores  python populate_mocks.py

srun --nodes=4 --tasks-per-node=1 --cpus-per-task=64 --cpu-bind=cores  python populate_mocks.py

/global/cfs/cdirs/des/dirac_sims/original_files/runsC/run001> run.00100.lightcone.npy
/pscratch/sd/m/mgatti/Dirac

'''

  






def make_maps(uuu):
 
    [seed,folder] = uuu
    # params files:
    #seed = 120
    # label
    if seed <10:
        mock_number = '00{0}'.format(seed)
    elif (seed>=10) & (seed<100):
        mock_number = '0{0}'.format(seed)
    elif (seed>=100):
        mock_number = '{0}'.format(seed)
    try: 
        if not os.path.exists(config['output']+'/runs{0}/'.format(folder,mock_number)):
            os.mkdir(config['output']+'/runs{0}/'.format(folder,mock_number))
    except:
        pass  

    try: 
        if not os.path.exists(config['output']+'/runs{0}/run{1}'.format(folder,mock_number)):
            os.mkdir(config['output']+'/runs{0}/run{1}'.format(folder,mock_number))
    except:
        pass  

    path_folder = config['path_mocks']+'/runs{0}/'.format(folder)+'/run{1}/'.format(folder,mock_number)
    path_folder_output = config['output']+'/runs{0}//run{1}/'.format(folder,mock_number)


  
    f = open(('params_run_1_Niall_{0}.txt'.format(folder)),'r')
    om_ = []
    h_ = []
    for i,f_ in enumerate(f):
        if i>0:
            om_.append(float(f_.split(',')[0]))
            h_.append(float(f_.split(',')[4]))

    om = om_[seed-1]
    h = h_[seed-1]*100.*u.km / u.s / u.Mpc

    # read redshift information ********************************************************************
 
    build_z_values_file(path_folder,'run',path_folder_output)

    resume = dict()
    resume['Step'] = []
    resume['z_far'] = []
    resume['z_near'] = []
    resume['delta_z'] = []
    resume['cmd_far'] = []
    resume['cmd_near'] = []
    resume['delta_cmd'] = []

    fil_ = open(path_folder_output+'/z_values.txt')
    for z__,z_ in enumerate(fil_):

            if z__>0:
                mute = np.array(z_.split(',')).astype(np.float)
                resume['Step'].append(mute[0])
                resume['z_far'].append(mute[1])
                resume['z_near'].append(mute[2])
                resume['delta_z'].append(mute[3])
                resume['cmd_far'].append(mute[4]/h_[seed-1])
                resume['cmd_near'].append(mute[5]/h_[seed-1])
                resume['delta_cmd'].append(mute[6]/h_[seed-1])

    run_param_sample = np.genfromtxt('params_run_1_Niall_{0}.txt'.format(folder),delimiter=',')[seed-1]
    w0 = run_param_sample[2]
  

    ###params = np.genfromtxt(path_folder +'transfer_function_cosmology.txt', dtype=None, skip_header=2)
###
    ###Omega_c=params[:,2][np.where(params[:,0] == b'Omega_cdm')].astype(float)[0]
    ###w0 = params[:,2][np.where(params[:,0] == b'w0_fld')].astype(float)[0]
    ###Omega_b=params[:,2][np.where(params[:,0] == b'Omega_b')].astype(float)[0]
    ####h=params[:,2][np.where(params[:,0] == b'h')].astype(float)[0]
    ###A_s=params[:,2][np.where(params[:,0] == b'A_s')].astype(float)[0]
    ###n_s= params[:,2][np.where(params[:,0] == b'n_s')].astype(float)[0]
    ###m_nu=0.06
    ###Omega_k=0.
    ###Omega_fld = params[:,2][np.where(params[:,0] == b'Omega_fld')].astype(float)[0]

   
    #  ********************************************************************************************

    # let's read raw particle numbers and make lens files:
    print ('save lenses')
    for s in frogress.bar(range(len(resume['Step']))):
        
        
        step_ = np.int(np.float(resume['Step'][s]))
      
        if step_ <10:
            zz = copy.copy('0000'+str(step_))
      
        elif (step_>=10) & (step_<100):
  
            zz =  copy.copy('000'+str(step_))
        elif (step_>=100):
                        
            zz =  copy.copy('00'+str(step_))

        path_ = path_folder_output+'/lens_{0}_{1}.fits'.format(np.int(step_),config['nside_intermediate'])
 
        if not os.path.exists(path_):
        
            if os.path.exists(path_folder+'/run.'+zz+'.lightcone.npy'.format(zz)):
                shell_ =np.load(path_folder+'/run.'+zz+'.lightcone.npy'.format(zz),allow_pickle=True)
                shell_ =  (shell_-np.mean(shell_))/np.mean(shell_)
                shell_ = hp.ud_grade(shell_,nside_out=config['nside_intermediate'])
                fits_f = Table()
                fits_f['T'] = shell_
                if os.path.exists(path_):
                    os.remove(path_)
                 
                fits_f.write(path_)

            
    
    # this is at reasonably high redshift...
    path_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(20,config['nside_intermediate'])



    kappa_pref_evaluated = brk.kappa_prefactor(h,om, length_unit = 'Mpc')
    
    
    z_near = np.array(resume['z_near'][::-1])
    z_far = np.array(resume['z_far'][::-1])
    z_bin_edges = np.hstack([z_near,z_far[-1]])
 
    z_bin_edges[0] = 1e-6
    from astropy.cosmology import FlatLambdaCDM,wCDM

    cosmology = wCDM(H0= h,
                 Om0=om,#mega_fld,
                 Ode0=1-om,#Omega_fld,
                 w0=w0)

    comoving_edges =  cosmology.comoving_distance(z_bin_edges)
    z_centre = np.array([z_at_value(cosmology.comoving_distance, 0.5*(comoving_edges[i]+comoving_edges[i+1]))  for i in range(len(comoving_edges)-1)])
    
    comoving_edges =  [cosmology.comoving_distance(x_) for x_ in np.array((z_bin_edges))]

    

    
    un_ = comoving_edges[:i][0].unit
    comoving_edges = np.array([c.value for c in comoving_edges])
    comoving_edges = comoving_edges*un_


    if not os.path.exists(path_):
        
    # let's make kappa files *************************************************************************



        overdensity_array = [np.zeros(hp.nside2npix(config['nside_intermediate']))]
        print ('load lens')
        for i in frogress.bar(range(len(resume['Step']))):
            try:
                shell = np.int(resume['Step'][-2])-i
                path_ = path_folder_output+'/lens_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
                m_ = pf.open(path_)
                overdensity_array.append(m_[1].data['T'])
            except:
                if i !=0:
                    overdensity_array.append(np.zeros(hp.nside2npix(config['nside_intermediate'])))


        print ('done ++++')
        overdensity_array = np.array(overdensity_array)
    

        from bornraytrace import lensing
        kappa_lensing = np.copy(overdensity_array)*0.

        print ('doing kappa')

      
        for i in frogress.bar(np.arange(1,kappa_lensing.shape[0]+1)):
                try:
                    kappa_lensing[i-1] = lensing.raytrace(h, om,
                                             overdensity_array=overdensity_array[:i].T,
                                             a_centre=1./(1.+z_bin_edges[1:(i+1)]),
                                             comoving_edges=comoving_edges[:(i+1)])
                except:
                    pass
    
        print ('done kappa')
        for i in range(kappa_lensing.shape[0]):
          
                shell = np.int(resume['Step'][-2])-i
                path_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
                if os.path.exists(path_):
                    os.remove(path_)
                fits_f = Table()
                fits_f['T'] = kappa_lensing[i]
                fits_f.write(path_)
            
        print ('save kappa done')
        
        del overdensity_array
        del kappa_lensing
        gc.collect()

        
        
    # lets make g1,g2 and IA files ********************
    print ('doing g1g2 ')
    for i in frogress.bar(range(len(resume['Step']))):
        shell = np.int(resume['Step'][-2])-i
        try:
            path_gg = path_folder_output+'/gg_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
            m = pf.open(path_gg)
            m = m[1].data['g1']
            

            
        except:
            if os.path.exists(path_gg):
                    os.remove(path_gg)
            try:
                
                path_ = path_folder_output+'/lens_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
                d = pf.open(path_)
                path_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
                k = pf.open(path_)
                path_ = path_folder_output+'/gg_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
          

                g1_, g2_ = gk_inv(k[1].data['T']-np.mean(k[1].data['T']),k[1].data['T']*0.,config['nside_intermediate'],config['nside_intermediate']*2)
                g1_IA, g2_IA = gk_inv(d[1].data['T']-np.mean(d[1].data['T']),k[1].data['T']*0.,config['nside_intermediate'],config['nside_intermediate']*2)

                fits_f = Table()
                fits_f['g1'] = g1_
                fits_f['g2'] = g2_
                fits_f['g1_IA'] = g1_IA
                fits_f['g2_IA'] = g2_IA
                fits_f.write(path_)
            except:
                pass

    
    # let's make integrated shear maps ***************************************************************************************************
    
    c1 = (5e-14 * (u.Mpc**3.)/(u.solMass * u.littleh**2) ) 
    c1_cgs = (c1* ((u.littleh/(cosmology.H0.value/100))**2.)).cgs
    rho_c1 = (c1_cgs*cosmology.critical_density(0)).value



    
    
    
    # load redshift distributions  **************
    
    import timeit
    st = timeit.default_timer()
    if docat:
        import timeit
        
        print ('LOADING catalog')
        
        sources_cat = dict()


        ## ROTATIONS ---------------------------------------------------------
        for rot in range(4):

            sources_cat[rot] = dict()
            mu = pf.open(config['2PT_FILE'])
            random_rel = np.random.randint(0,6000,1)[0]
            redshift_distributions_sources = {'z':None,'bins':dict()}
            redshift_distributions_sources['z'] = mu[8+random_rel].data['Z_MID']
            for ix in config['sources_bins']:
                redshift_distributions_sources['bins'][ix] = mu[8+random_rel].data['BIN{0}'.format(ix)]
            mu = None


            config['A_IA'] = np.random.randint(-300000,300000,1)[0]/100000
            config['eta_IA'] = np.random.randint(-500000,500000,1)[0]/100000
            config['A_IA'] = np.random.randint(-300000,300000,1)[0]/100000
            #config['eta_IA'] = 0.
            config['z0_IA'] = 0.67
            BIAS_SC = 1.

            g1_tomo = dict()
            g2_tomo = dict()
            d_tomo = dict()
            nz_kernel_sample_dict = dict()

            print ('LOADING, {0}'.format(rot))
            for tomo_bin in config['sources_bins']:
                g1_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside_intermediate']))
                g2_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside_intermediate']))
                d_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside_intermediate']))
                redshift_distributions_sources['bins'][tomo_bin][250:] = 0.
                nz_sample = brk.recentre_nz(np.array(z_bin_edges).astype('float'),  redshift_distributions_sources['z'],  redshift_distributions_sources['bins'][tomo_bin] )
                nz_kernel_sample_dict[tomo_bin] = nz_sample*(z_bin_edges[1:]-z_bin_edges[:-1])
               

                for i in frogress.bar(range(2,len(comoving_edges)-1)):

                    try:
                        shell = np.int(resume['Step'][-2])-i
                        path_ = path_folder_output+'/lens_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
                        pathk_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
                        pathgg_ =  path_folder_output+'/gg_{0}_{1}.fits'.format(shell,config['nside_intermediate'])




                        k_ = pf.open(pathgg_)
                        k_real_ = pf.open(pathk_)
                        d_ = pf.open(path_)
                        IA_f = iaa.F_nla(z_centre[i], cosmology.Om0, rho_c1=rho_c1,A_ia = config['A_IA'], eta=config['eta_IA'], z0=config['z0_IA'],  lbar=0., l0=1e-9, beta=0.)
                        #print ((k_[1].data['T']))



                        g1_tomo[tomo_bin]  +=  ((1.+BIAS_SC*(d_[1].data['T']))*(k_[1].data['g1']+k_[1].data['g1_IA']*IA_f))*nz_kernel_sample_dict[tomo_bin][i]
                        g2_tomo[tomo_bin]  +=  ((1.+BIAS_SC*(d_[1].data['T']))*(k_[1].data['g2']+k_[1].data['g2_IA']*IA_f))*nz_kernel_sample_dict[tomo_bin][i]
                        d_tomo[tomo_bin] +=  (1.+BIAS_SC*d_[1].data['T'])*nz_kernel_sample_dict[tomo_bin][i]
                    except:
                        pass
                
                if rot ==0:
                    g1_tomo[tomo_bin] = hp.ud_grade( copy.copy(g1_tomo[tomo_bin]),nside_out=config['nside2'])
                    g2_tomo[tomo_bin] = hp.ud_grade( copy.copy(g2_tomo[tomo_bin]),nside_out=config['nside2'])
                    d_tomo[tomo_bin] =  hp.ud_grade( copy.copy(d_tomo[tomo_bin]),nside_out=config['nside2'])
                elif rot ==1:

                    g1_tomo[tomo_bin] = rotate_map_approx(hp.ud_grade(g1_tomo[tomo_bin],nside_out=config['nside2']),[ 180 ,0 , 0], flip=False,nside = config['nside2'])
                    g2_tomo[tomo_bin] = rotate_map_approx(hp.ud_grade(g2_tomo[tomo_bin],nside_out=config['nside2']),[ 180 ,0 , 0], flip=False,nside = config['nside2'])
                    d_tomo[tomo_bin] =  rotate_map_approx(hp.ud_grade(d_tomo[tomo_bin] ,nside_out=config['nside2']),[ 180 ,0 , 0], flip=False ,nside = config['nside2'])
                elif rot ==2:

                    g1_tomo[tomo_bin] = rotate_map_approx(hp.ud_grade(g1_tomo[tomo_bin],nside_out=config['nside2']),[ 90 ,0 , 0], flip=True,nside = config['nside2'] )
                    g2_tomo[tomo_bin] = rotate_map_approx(hp.ud_grade(g2_tomo[tomo_bin],nside_out=config['nside2']),[ 90 ,0 , 0], flip=True,nside = config['nside2'] )
                    d_tomo[tomo_bin] =  rotate_map_approx(hp.ud_grade(d_tomo[tomo_bin] ,nside_out=config['nside2']),[ 90 ,0 , 0], flip=True ,nside = config['nside2'] )
                elif rot ==3:

                    g1_tomo[tomo_bin] = rotate_map_approx(hp.ud_grade(g1_tomo[tomo_bin],nside_out=config['nside2']),[ 270 ,0 , 0], flip=True,nside = config['nside2'] )
                    g2_tomo[tomo_bin] = rotate_map_approx(hp.ud_grade(g2_tomo[tomo_bin],nside_out=config['nside2']),[ 270 ,0 , 0], flip=True,nside = config['nside2'] )
                    d_tomo[tomo_bin] =  rotate_map_approx(hp.ud_grade(d_tomo[tomo_bin] ,nside_out=config['nside2']),[ 270 ,0 , 0], flip=True ,nside = config['nside2'] )

                    
            print ('done loading')      
                    


            nuis = dict()
            nuis['m'] = np.zeros(4) 
            for tomo_bin in config['sources_bins']:
                m_ = 1+config['m_sources'][tomo_bin-1]
                s_ = config['ms_sources'][tomo_bin-1]
                m_1 = np.random.normal(m_,s_,1)[0]
                nuis['m'][tomo_bin-1] = m_1









            config['nside2'] = 512 
            depth_weigth = np.load('/global/cfs/cdirs/des/mass_maps/Maps_final/depth_maps_Y3_{0}_numbdensity.npy'.format(config['nside2']),allow_pickle=True).item()

            for tomo_bin in config['sources_bins']:

                sources_cat[rot][tomo_bin] = dict()

                
                mcal_catalog = load_obj('/global/cfs/cdirs/des/mass_maps/Maps_final/data_catalogs_weighted_{0}'.format(tomo_bin-1))


                pix_ = convert_to_pix_coord(mcal_catalog['ra'], mcal_catalog['dec'], nside=config['nside2'])
                mask = np.in1d(np.arange(hp.nside2npix(config['nside2'])),pix_)


                # generate ellipticities ***********************************
                df2 = pd.DataFrame(data = {'w':mcal_catalog['w'] ,'pix_':pix_},index = pix_)
                #df2 = df2.groupby(df2.columns.tolist(),as_index=False).size()
                
                #print (len(df2))
                # poisson sample the weight map ------
                nn = np.random.poisson(depth_weigth[tomo_bin-1])

                nn[~mask]= 0


                count = 0

                nnmaxx = max(nn)
                    
                for count in range(nnmaxx):

                    if count %2 ==0:
                        df3 = df2.sample(frac=1)
                        df4 = df3.drop_duplicates('pix_',keep ='first').sort_index()
                    else:
                        df4 = df3.drop_duplicates('pix_',keep ='last').sort_index()
                        
                    pix_valid = np.arange(len(nn))[nn>0]
                    df3 = df4.loc[np.unique(pix_valid)]
                    if count == 0:
                        w = df3['w']
                        pix = df3['pix_']
                    else:
                        w = np.hstack([w,df3['w']])
                        pix = np.hstack([pix,df3['pix_']]) 
                    nn -= 1
       
                del df2
                del df3
                gc.collect()
                e1,e2 = random_draw_ell_from_w(w,mcal_catalog['w'],mcal_catalog['e1'],mcal_catalog['e2'])



                del mcal_catalog
                gc.collect()


                f = 1./np.sqrt(d_tomo[tomo_bin]/np.sum(nz_kernel_sample_dict[tomo_bin]))



             
                f = f[pix]


                # ++++++++++++++++++++++
 
                n_map_sc = np.zeros(hp.nside2npix(config['nside2']))

                unique_pix, idx, idx_rep = np.unique(pix, return_index=True, return_inverse=True)

       
                n_map_sc[unique_pix] += np.bincount(idx_rep, weights=w/f**2)

                g1_ = g1_tomo[tomo_bin][pix]
                g2_ = g2_tomo[tomo_bin][pix]


                es1,es2 = apply_random_rotation(e1/f, e2/f)
                es1a,es2a = apply_random_rotation(e1/f, e2/f)


                x1_sc,x2_sc = addSourceEllipticity({'shear1':g1_,'shear2':g2_},{'e1':es1,'e2':es2},es_colnames=("e1","e2"))


                e1r_map = np.zeros(hp.nside2npix(config['nside2']))
                e2r_map = np.zeros(hp.nside2npix(config['nside2']))

                e1r_map0 = np.zeros(hp.nside2npix(config['nside2']))
                e2r_map0 = np.zeros(hp.nside2npix(config['nside2']))

                g1_map = np.zeros(hp.nside2npix(config['nside2']))
                g2_map = np.zeros(hp.nside2npix(config['nside2']))

                unique_pix, idx, idx_rep = np.unique(pix, return_index=True, return_inverse=True)




                e1r_map[unique_pix] += np.bincount(idx_rep, weights=es1*w)
                e2r_map[unique_pix] += np.bincount(idx_rep, weights=es2*w)

                e1r_map0[unique_pix] += np.bincount(idx_rep, weights=es1a*w)
                e2r_map0[unique_pix] += np.bincount(idx_rep, weights=es2a*w)


                g1_map[unique_pix] += np.bincount(idx_rep, weights= g1_*w)
                g2_map[unique_pix] += np.bincount(idx_rep, weights= g2_*w)







                mask_sims = n_map_sc != 0.
                e1r_map[mask_sims]  = e1r_map[mask_sims]/(n_map_sc[mask_sims])
                e2r_map[mask_sims] =  e2r_map[mask_sims]/(n_map_sc[mask_sims])
                e1r_map0[mask_sims]  = e1r_map0[mask_sims]/(n_map_sc[mask_sims])
                e2r_map0[mask_sims] =  e2r_map0[mask_sims]/(n_map_sc[mask_sims])
                g1_map[mask_sims]  = g1_map[mask_sims]/(n_map_sc[mask_sims])
                g2_map[mask_sims] =  g2_map[mask_sims]/(n_map_sc[mask_sims])




                #EE,BB,_   =  g2k_sphere((g1_map+e1r_map0)*nuis['m'][tomo_bin-1], (g2_map+e2r_map0)*nuis['m'][tomo_bin-1], mask_sims, nside=config['nside2'], lmax=config['nside2']*2 ,nosh=True)
               # EEn,BBn,_ =  g2k_sphere(e1r_map*nuis['m'][tomo_bin-1], e2r_map*nuis['m'][tomo_bin-1], mask_sims, nside=config['nside2'], lmax=config['nside2']*2 ,nosh=True)
                #sources_cat[rot][tomo_bin] = {'kE':EE,'kE_noise':EEn,'mask':mask_sims}
                
                
                e1_ = ((g1_map+e1r_map0)*nuis['m'][tomo_bin-1])[mask_sims]
                e2_ = ((g2_map+e2r_map0)*nuis['m'][tomo_bin-1])[mask_sims]
                e1n_ = ( e1r_map*nuis['m'][tomo_bin-1])[mask_sims]
                e2n_ = ( e2r_map*nuis['m'][tomo_bin-1])[mask_sims]
                idx_ = np.arange(len(mask_sims))[mask_sims]
                
                sources_cat[rot][tomo_bin] = {'e1':e1_,'e2':e2_,'e1n':e1n_,'e2n':e2n_,'pix':idx_}

                #sources_cat[tomo_bin]['k_orig'] = EE_orig
               # del mcal_catalog[tomo_bin-1]
                
            nuis['hyperrank_rel'] = random_rel
            nuis['A_IA'] = config['A_IA'] 
            nuis['E_IA'] = config['eta_IA'] 
            sources_cat[rot]['nuisance_params'] = nuis
            sources_cat[rot]['config'] = copy.deepcopy(config)






        np.save(config['output']+'/shear_maps_runs{0}_run{1}_nside{2}_noiserel{3}'.format(folder,mock_number,config['nside2'],config['noise_rel']),sources_cat)
        del sources_cat
        gc.collect()

  
    






if __name__ == '__main__':
    
    
    docat=True
    folders = ['C','E','I','J','K','L','M','N','O','P','Q','R','S']
    #folders = ['K','L','M','N','O','P','Q','R','S']
    #folders = ['R','S']
    
    # 10 C,4 I,  20M 32o 32p 20Q
    #folders = ['P']#,'R','S']
    runstodo=[]
    for folder_ in folders:
        config = dict()
        config['noise_rel'] = 14
        config['2PT_FILE'] = '//global/cfs/cdirs//des/www/y3_chains/data_vectors/2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate_6000HR.fits'
        config['nside_intermediate'] = 512
        config['nside'] = 512
        config['path_mocks'] = '/global/cfs/cdirs/des/dirac_sims/original_files/'
        #/global/cfs/cdirs/des/mgatti/Dirac
        config['output'] = '/pscratch/sd/m/mgatti/Dirac/' #/global/cfs/cdirs/des/dirac_sims/derived_products/'
        config['sources_bins'] = [1,2,3,4]#,2,3,4]#,2,3,4] #1,2,3,4

        config['m_sources'] = [-0.002,-0.017,-0.029,-0.038]#,[-0.006,-.01,-0.026,-0.032]
        config['ms_sources'] = [0.0091,0.0078,0.0076,0.0076]


        #make folder:
        try:
            if not os.path.exists(config['output']+'/runs{0}/'.format(folder_)):
                os.mkdir(config['output']+'/runs{0}/'.format(folder_))
        except:
            pass

        # figure out how many realisations in the folder **************************************************
        import numpy as np
        import glob
        files = glob.glob(config['path_mocks']+'/runs{0}/'.format(folder_)+'/*')
        rel =[]
        for file in files:
            try:
                rel.append(file.split('run')[2].strip('.tar.gz') )
            except:
                pass
    
        #rel = [file.split('run')[2].strip('.tar.gz') for file in files]
        rel_ = []
        for r in rel:
            try:
                rel_.append(float(r))
            except:
                pass
        config['n_mocks'] = len(np.unique(rel_))
        #**************************************************************************************************
        config['nside2'] = 512
        


        for seed in range(1,config['n_mocks']+1):
            if seed <10:
                mock_number = '00{0}'.format(seed)
            elif (seed>=10) & (seed<100):
                mock_number = '0{0}'.format(seed)
            elif (seed>=100):
                mock_number = '{0}'.format(seed)
            path = config['output']+'/shear_maps_runs{0}_run{1}_nside{2}_noiserel{3}'.format(folder_,mock_number,config['nside2'],config['noise_rel'])+'.npy'
            #path = '/global/cscratch1/sd/ucapnje/DES_DiRAC/runs_kappa/runs{0}/Ashear_maps_'.format(folder)+str(mock_number)+'.npy'
            if not os.path.exists(path):
                runstodo.append([seed,folder_])
    run_count=0
    print (len(runstodo))
        
        
    from mpi4py import MPI 
### 
    while run_count<len(runstodo):
        comm = MPI.COMM_WORLD
        print("Hello! I'm rank %d from %d running in total..." % (comm.rank, comm.size))
        if (run_count+comm.rank)<len(runstodo):
            try:
                make_maps(runstodo[run_count+comm.rank])
            except:
                print ('failed ',runstodo[run_count+comm.rank])
             #   pass
        run_count+=comm.size
        comm.bcast(run_count,root = 0)
        comm.Barrier() 
       
#   from mpi4py import MPI 
 
  #while run_count<len(runstodo):
#
  #    #print("Hello! I'm rank %d from %d running in total..." % (comm.rank, comm.size))
  #    if (run_count)<len(runstodo):
  #        #try:
  #        make_maps(runstodo[run_count])
  #        #except:
  #        #    print ('failed ',runstodo[run_count+comm.rank])
  #        #    pass
  #    run_count+=1 #comm.size
  #   # comm.bcast(run_count,root = 0)
  #    #comm.Barrier() 
  ##
