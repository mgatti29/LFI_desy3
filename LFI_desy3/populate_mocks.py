import bornraytrace
from bornraytrace import lensing
from bornraytrace import lensing as brk
from bornraytrace import intrinsic_alignments as iaa
import gc
import pyfits as pf
import numpy as np
import healpy as hp
import pandas as pd
import pickle
import os
import copy
import scipy
from scipy.interpolate import interp1d
import timeit
from astropy import units as u
from astropy.table import Table  
from astropy.cosmology import z_at_value
from Moments_analysis import convert_to_pix_coord, IndexToDeclRa, apply_random_rotation, addSourceEllipticity, gk_inv
import frogress
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
srun --nodes=4 --tasks-per-node=64 python populate_mocks_corr.py


'''

  






def make_maps(input_):
 
    # reads seed of the mock and folder
    [seed,folder] = input_

    #convert to the right format the seed
    if seed <10:
        mock_number = '00{0}'.format(seed)
    elif (seed>=10) & (seed<100):
        mock_number = '0{0}'.format(seed)
    elif (seed>=100):
        mock_number = '{0}'.format(seed)
        
    #makes folder if it doesn't exist
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

    # path to folders
    path_folder = config['path_mocks']+'/runs{0}/'.format(folder)+'/run{1}/'.format(folder,mock_number)
    path_folder_output = config['output']+'/runs{0}//run{1}/'.format(folder,mock_number)


    
    # this reads the cosmological parameter of the simulations
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
  

   
    #  ********************************************************************************************

    # let's read raw particle numbers and make lens files:
    path_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(30,config['nside'])
    if not os.path.exists(path_):
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
                    
            path_ = path_folder_output+'/lens_{0}_{1}.fits'.format(np.int(step_),config['nside'])

            if not os.path.exists(path_):

                if os.path.exists(path_folder+'/run.'+zz+'.lightcone.npy'.format(zz)):
                    shell_ =np.load(path_folder+'/run.'+zz+'.lightcone.npy'.format(zz),allow_pickle=True)
                    shell_ =  (shell_-np.mean(shell_))/np.mean(shell_)
                    shell_ = hp.ud_grade(shell_,nside_out=config['nside'])
                    fits_f = Table()
                    fits_f['T'] = shell_
                    if os.path.exists(path_):
                        os.remove(path_)

                    fits_f.write(path_)

            
    
    # this is at reasonably high redshift...
    path_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(30,config['nside'])


    
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
                path_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(shell,config['nside'])
                if os.path.exists(path_):
                    os.remove(path_)
                fits_f = Table()
                fits_f['T'] = hp.ud_grade(kappa_lensing[i],nside_out = config['nside'])
                fits_f.write(path_)
            
            
                path_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
                if os.path.exists(path_):
                    os.remove(path_)
                fits_f = Table()
                fits_f['T'] = hp.ud_grade(kappa_lensing[i],nside_out =config['nside_intermediate'])
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
            path_gg = path_folder_output+'/gg_{0}_{1}.fits'.format(shell,config['nside'])
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


                g1_, g2_ = gk_inv(k[1].data['T']-np.mean(k[1].data['T']),k[1].data['T']*0.,config['nside_intermediate'],config['nside_intermediate']*2)
                g1_IA, g2_IA = gk_inv(d[1].data['T']-np.mean(d[1].data['T']),k[1].data['T']*0.,config['nside_intermediate'],config['nside_intermediate']*2)

                path_ = path_folder_output+'/gg_{0}_{1}.fits'.format(shell,config['nside'])
          
                fits_f = Table()
                fits_f['g1'] = hp.ud_grade( g1_,nside_out =config['nside'])
                fits_f['g2'] = hp.ud_grade(g2_,nside_out =config['nside'])
                fits_f['g1_IA'] = hp.ud_grade(g1_IA,nside_out =config['nside'])
                fits_f['g2_IA'] = hp.ud_grade(g2_IA,nside_out =config['nside'])
                fits_f.write(path_)
            except:
                pass

            path_ = path_folder_output+'/lens_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
            os.system('rm {0}'.format(path_))
            path_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(shell,config['nside_intermediate'])
            os.system('rm {0}'.format(path_))
          
    
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
            if not SC:
                BIAS_SC = 0.
            else:
                BIAS_SC = 1.

            g1_tomo = dict()
            g2_tomo = dict()
            d_tomo = dict()
            nz_kernel_sample_dict = dict()

            print ('LOADING, {0}'.format(rot))
            for tomo_bin in config['sources_bins']:
                g1_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
                g2_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
                d_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
                redshift_distributions_sources['bins'][tomo_bin][250:] = 0.
                nz_sample = brk.recentre_nz(np.array(z_bin_edges).astype('float'),  redshift_distributions_sources['z'],  redshift_distributions_sources['bins'][tomo_bin] )
                nz_kernel_sample_dict[tomo_bin] = nz_sample*(z_bin_edges[1:]-z_bin_edges[:-1])
               

                for i in frogress.bar(range(2,len(comoving_edges)-1)):

                    try:
                        shell = np.int(resume['Step'][-2])-i
                        path_ = path_folder_output+'/lens_{0}_{1}.fits'.format(shell,config['nside'])
                        pathk_ = path_folder_output+'/kappa_{0}_{1}.fits'.format(shell,config['nside'])
                        pathgg_ =  path_folder_output+'/gg_{0}_{1}.fits'.format(shell,config['nside'])




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
                
                #if rot ==0:
                g1_tomo[tomo_bin] = hp.ud_grade( copy.copy(g1_tomo[tomo_bin]),nside_out=config['nside2'])
                g2_tomo[tomo_bin] = hp.ud_grade( copy.copy(g2_tomo[tomo_bin]),nside_out=config['nside2'])
                d_tomo[tomo_bin] =  hp.ud_grade( copy.copy(d_tomo[tomo_bin]),nside_out=config['nside2'])
                
                
               


            print ('done loading')      
                    


            nuis = dict()
            nuis['m'] = np.zeros(4) 
            for tomo_bin in config['sources_bins']:
                m_ = 1+config['m_sources'][tomo_bin-1]
                s_ = config['ms_sources'][tomo_bin-1]
                m_1 = np.random.normal(m_,s_,1)[0]
                nuis['m'][tomo_bin-1] = m_1







            # let's make DES mock catalogs

            config['nside2'] = 512 

            for tomo_bin in config['sources_bins']:

                sources_cat[rot][tomo_bin] = dict()
                
                # load into memory the des y3 mock catalogue
                mcal_catalog = load_obj('/global/cfs/cdirs/des/mass_maps/Maps_final/data_catalogs_weighted_{0}'.format(tomo_bin-1))
                
                pix_ = convert_to_pix_coord(mcal_catalog['ra'], mcal_catalog['dec'], nside=config['nside2'])
                
                e1 = mcal_catalog['e1']
                e2 = mcal_catalog['e2']
                w  = mcal_catalog['w']
                
                # we can cut 4 des y3 footprint out of the full sky. this is controlled by the 'rot' parameter
                if rot ==0:
                    rot_angles = [0, 0, 0]
                    flip=False
                    rotu = hp.rotator.Rotator(rot=rot_angles, deg=True)
                    alpha, delta = hp.pix2ang(config['nside2'],pix_)
                    rot_alpha, rot_delta = rotu(alpha, delta)
                    if not flip:
                        pix = hp.ang2pix(config['nside2'], rot_alpha, rot_delta)
                    else:
                        pix = hp.ang2pix(config['nside2'], np.pi-rot_alpha, rot_delta)
             
                if rot ==1:
                    rot_angles = [180, 0, 0]
                    flip=False
                    rotu = hp.rotator.Rotator(rot=rot_angles, deg=True)
                    alpha, delta = hp.pix2ang(config['nside2'],pix_)
                    rot_alpha, rot_delta = rotu(alpha, delta)
                    if not flip:
                        pix = hp.ang2pix(config['nside2'], rot_alpha, rot_delta)
                    else:
                        pix = hp.ang2pix(config['nside2'], np.pi-rot_alpha, rot_delta)
              
                if rot ==2:
                    rot_angles = [90, 0, 0]
                    flip=True
                    rotu = hp.rotator.Rotator(rot=rot_angles, deg=True)
                    alpha, delta = hp.pix2ang(config['nside2'],pix_)
                    rot_alpha, rot_delta = rotu(alpha, delta)
                    if not flip:
                        pix = hp.ang2pix(config['nside2'], rot_alpha, rot_delta)
                    else:
                        pix = hp.ang2pix(config['nside2'], np.pi-rot_alpha, rot_delta)
                
                if rot ==3:
                    rot_angles = [270, 0, 0]
                    flip=True
                    rotu = hp.rotator.Rotator(rot=rot_angles, deg=True)
                    alpha, delta = hp.pix2ang(config['nside2'],pix_)
                    rot_alpha, rot_delta = rotu(alpha, delta)
                    if not flip:
                        pix = hp.ang2pix(config['nside2'], rot_alpha, rot_delta)
                    else:
                        pix = hp.ang2pix(config['nside2'], np.pi-rot_alpha, rot_delta)
                      

                del mcal_catalog
                gc.collect() 
                
                # the factor f controls the source clustering effect on shape noise
                f = 1./np.sqrt(d_tomo[tomo_bin]/np.sum(nz_kernel_sample_dict[tomo_bin]))

                f = f[pix]

                if not SC:
                    f = 1.

                # ++++++++++++++++++++++
 
                n_map = np.zeros(hp.nside2npix(config['nside2']))
                n_map_sc = np.zeros(hp.nside2npix(config['nside2']))

                unique_pix, idx, idx_rep = np.unique(pix, return_index=True, return_inverse=True)

       
                n_map_sc[unique_pix] += np.bincount(idx_rep, weights=w/f**2)
                n_map[unique_pix] += np.bincount(idx_rep, weights=w)

                g1_ = g1_tomo[tomo_bin][pix]
                g2_ = g2_tomo[tomo_bin][pix]


                es1,es2 = apply_random_rotation(e1/f, e2/f)
                es1_ref,es2_ref = apply_random_rotation(e1, e2)
                es1a,es2a = apply_random_rotation(e1/f, e2/f)


                x1_sc,x2_sc = addSourceEllipticity({'shear1':g1_,'shear2':g2_},{'e1':es1,'e2':es2},es_colnames=("e1","e2"))


                e1r_map = np.zeros(hp.nside2npix(config['nside2']))
                e2r_map = np.zeros(hp.nside2npix(config['nside2']))

                e1r_map0 = np.zeros(hp.nside2npix(config['nside2']))
                e2r_map0 = np.zeros(hp.nside2npix(config['nside2']))

                e1r_map0_ref = np.zeros(hp.nside2npix(config['nside2']))
                e2r_map0_ref = np.zeros(hp.nside2npix(config['nside2']))

                g1_map = np.zeros(hp.nside2npix(config['nside2']))
                g2_map = np.zeros(hp.nside2npix(config['nside2']))

                unique_pix, idx, idx_rep = np.unique(pix, return_index=True, return_inverse=True)


                e1r_map[unique_pix] += np.bincount(idx_rep, weights=es1*w)
                e2r_map[unique_pix] += np.bincount(idx_rep, weights=es2*w)

                e1r_map0[unique_pix] += np.bincount(idx_rep, weights=es1a*w)
                e2r_map0[unique_pix] += np.bincount(idx_rep, weights=es2a*w)

                e1r_map0_ref[unique_pix] += np.bincount(idx_rep, weights=es1_ref*w)
                e2r_map0_ref[unique_pix] += np.bincount(idx_rep, weights=es2_ref*w)

                
                mask_sims = n_map_sc != 0.
                e1r_map[mask_sims]  = e1r_map[mask_sims]/(n_map_sc[mask_sims])
                e2r_map[mask_sims] =  e2r_map[mask_sims]/(n_map_sc[mask_sims])
                e1r_map0[mask_sims]  = e1r_map0[mask_sims]/(n_map_sc[mask_sims])
                e2r_map0[mask_sims] =  e2r_map0[mask_sims]/(n_map_sc[mask_sims])
                e1r_map0_ref[mask_sims]  = e1r_map0_ref[mask_sims]/(n_map[mask_sims])
                e2r_map0_ref[mask_sims] =  e2r_map0_ref[mask_sims]/(n_map[mask_sims])
                
                   
                
                var_ =  e1r_map0_ref**2+e2r_map0_ref**2

    
                #'''
                e1r_map[unique_pix]  *= 1/(np.sqrt(0.995*corr[tomo_bin-1])) * np.sqrt((1-coeff_kurtosis[tomo_bin-1]*var_))
                e2r_map[unique_pix]  *= 1/(np.sqrt(0.995*corr[tomo_bin-1])) * np.sqrt((1-coeff_kurtosis[tomo_bin-1]*var_))
                e1r_map0[unique_pix] *= 1/(np.sqrt(0.995*corr[tomo_bin-1])) * np.sqrt((1-coeff_kurtosis[tomo_bin-1]*var_))
                e2r_map0[unique_pix] *= 1/(np.sqrt(0.995*corr[tomo_bin-1])) * np.sqrt((1-coeff_kurtosis[tomo_bin-1]*var_))

                
                #'''
                g1_map[unique_pix] += np.bincount(idx_rep, weights= g1_*w)
                g2_map[unique_pix] += np.bincount(idx_rep, weights= g2_*w)








                g1_map[mask_sims]  = g1_map[mask_sims]/(n_map_sc[mask_sims])
                g2_map[mask_sims] =  g2_map[mask_sims]/(n_map_sc[mask_sims])
                
                e1_ = ((g1_map*nuis['m'][tomo_bin-1]+e1r_map0))[mask_sims]
                e2_ = ((g2_map*nuis['m'][tomo_bin-1]+e2r_map0))[mask_sims]
                e1n_ = ( e1r_map)[mask_sims]
                e2n_ = ( e2r_map)[mask_sims]
                idx_ = np.arange(len(mask_sims))[mask_sims]
                
                sources_cat[rot][tomo_bin] = {'e1':e1_,'e2':e2_,'e1n':e1n_,'e2n':e2n_,'pix':idx_}
   
            nuis['hyperrank_rel'] = random_rel
            nuis['A_IA'] = config['A_IA'] 
            nuis['E_IA'] = config['eta_IA'] 
            sources_cat[rot]['nuisance_params'] = nuis
            sources_cat[rot]['config'] = copy.deepcopy(config)






        np.save(config['output']+'/shear_maps_runs{0}_run{1}_nside{2}_noiserel{3}'.format(folder,mock_number,config['nside2'],config['noise_rel']),sources_cat)
        del sources_cat
        gc.collect()


    

corr = [1.0608,1.0295,1.0188,1.0115]
coeff_kurtosis = [0.1,0.05,0.036,0.036]

SC = True

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
        config['noise_rel'] = 9
        config['2PT_FILE'] = '//global/cfs/cdirs//des/www/y3_chains/data_vectors/2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate_6000HR.fits'
        config['nside_intermediate'] = 1024
        config['nside'] = 512
        config['path_mocks'] = '/global/cfs/cdirs/des/dirac_sims/original_files/'
        #/global/cfs/cdirs/des/mgatti/Dirac
        config['output'] = '/global/cfs/cdirs/des/mgatti/Dirac_mocks/' #/global/cfs/cdirs/des/dirac_sims/derived_products/'
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
    #make_maps(runstodo[0]) 
      
    

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
 
