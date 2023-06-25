import yaml
import pickle
import healpy as hp
import numpy as np
import os
from astropy.table import Table
import gc
import pyfits as pf
from Moments_analysis import g2k_sphere
import timeit
import os
import copy
from bornraytrace import lensing as brk
import numpy as np
from bornraytrace import intrinsic_alignments as iaa
import bornraytrace
from astropy.table import Table
import healpy as hp
import frogress
import pyfits as pf
from astropy.cosmology import z_at_value
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import cosmolopy.distance as cd
from scipy.interpolate import interp1d
import gc
import pandas as pd
import pickle
import multiprocessing
from functools import partial
from astropy.cosmology import FlatLambdaCDM,wCDM

def apply_random_rotation(e1_in, e2_in):
    np.random.seed() # CRITICAL in multiple processes !
    rot_angle = np.random.rand(len(e1_in))*2*np.pi #no need for 2?
    cos = np.cos(rot_angle)
    sin = np.sin(rot_angle)
    e1_out = + e1_in * cos + e2_in * sin
    e2_out = - e1_in * sin + e2_in * cos
    return e1_out, e2_out

def IndexToDeclRa(index, nside,nest= False):
    theta,phi=hp.pixelfunc.pix2ang(nside ,index,nest=nest)
    return -np.degrees(theta-np.pi/2.),np.degrees(phi)

def convert_to_pix_coord(ra, dec, nside=1024):
    """
    Converts RA,DEC to hpix coordinates
    """

    theta = (90.0 - dec) * np.pi / 180.
    phi = ra * np.pi / 180.
    pix = hp.ang2pix(nside, theta, phi, nest=False)

    return pix

def generate_randoms_radec(minra, maxra, mindec, maxdec, Ngen, raoffset=0):
    r = 1.0
    # this z is not redshift!
    zmin = r * np.sin(np.pi * mindec / 180.)
    zmax = r * np.sin(np.pi * maxdec / 180.)
    # parity transform from usual, but let's not worry about that
    phimin = np.pi / 180. * (minra - 180 + raoffset)
    phimax = np.pi / 180. * (maxra - 180 + raoffset)
    # generate ra and dec
    z_coord = np.random.uniform(zmin, zmax, Ngen)  # not redshift!
    phi = np.random.uniform(phimin, phimax, Ngen)
    dec_rad = np.arcsin(z_coord / r)
    # convert to ra and dec
    ra = phi * 180 / np.pi + 180 - raoffset
    dec = dec_rad * 180 / np.pi
    return ra, dec


def addSourceEllipticity(self,es,es_colnames=("e1","e2"),rs_correction=True,inplace=False):

		"""

		:param es: array of intrinsic ellipticities,

		"""

		#Safety check
		assert len(self)==len(es)

		#Compute complex source ellipticity, shear
		es_c = np.array(es[es_colnames[0]]+es[es_colnames[1]]*1j)
		g = np.array(self["shear1"] + self["shear2"]*1j)

		#Shear the intrinsic ellipticity
		e = es_c + g
		if rs_correction:
			e /= (1 + g.conjugate()*es_c)

		#Return
		if inplace:
			self["shear1"] = e.real
			self["shear2"] = e.imag
		else:
			return (e.real,e.imag)

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

def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, protocol=2)
        f.close()

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        mute =  pickle.load(f)
        f.close()
    return mute



def gk_inv(K,KB,nside,lmax):

    alms = hp.map2alm(K, lmax=lmax, pol=False)  # Spin transform!

    ell, emm = hp.Alm.getlm(lmax=lmax)

    kalmsE = alms/( 1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5)

    kalmsE[ell == 0] = 0.0


    alms = hp.map2alm(KB, lmax=lmax, pol=False)  # Spin transform!

    ell, emm = hp.Alm.getlm(lmax=lmax)

    kalmsB = alms/( 1. * ((ell * (ell + 1.)) / ((ell + 2.) * (ell - 1))) ** 0.5)

    kalmsB[ell == 0] = 0.0

    _,e1t,e2t = hp.alm2map([kalmsE,kalmsE,kalmsB] , nside=nside, lmax=lmax, pol=True)
    return e1t,e2t# ,r



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



def rotate_map_approx(mask, rot_angles, flip=False,nside = 2048):
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





def make_maps(seed):
    st = timeit.default_timer()

    # READ IN PARAMETERS ********************************************

    p,params_dict = seed

    # SET COSMOLOGY ************************************************
    config = dict()
    config['Om'] = params_dict['Omegam']
    config['sigma8'] =  params_dict['s8']
    config['ns'] =params_dict['ns']
    config['Ob'] = Ob
    config['h100'] = params_dict['h']

    config['nside_out'] = 512
    config['nside_intermediate'] = 1024
    config['nside'] = 512
    config['sources_bins'] = [1,2,3,4]
    config['dz_sources'] = [params_dict['dz1'],params_dict['dz2'],params_dict['dz3'],params_dict['dz4']]
    config['m_sources'] = [params_dict['m1'],params_dict['m2'],params_dict['m3'],params_dict['m4']]
    config['A_IA'] = params_dict['A']
    config['eta_IA'] = params_dict['E']
    config['f'] = params_dict['f']
    config['z0_IA'] = 0.67
    config['2PT_FILE'] = '/global/homes/m/mgatti/Mass_Mapping/HOD/PKDGRAV_CODE//2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate_6000HR.fits'

    rot = params_dict['rot']

    cosmo1 = {'omega_M_0': config['Om'],
     'omega_lambda_0':1-config['Om'],
     'omega_k_0':0.0,
     'omega_b_0' : config['Ob'],
     'h':config['h100'],
     'sigma_8' : config['sigma8'],
     'n': config['ns']}

    print ('h100',config['h100'])
    cosmology = wCDM(H0= config['h100']*u.km / u.s / u.Mpc,
                 Om0=config['Om'],#mega_fld,
                 Ode0=1-config['Om'],#Omega_fld,
                 w0=w0)



    '''
    Now the code will try to read in particle counts and make
    kappa,e1,e2 maps. If they are already there, this phase will be skipped
    '''

    # path where kappa/g1/g2 maps are stored
    base = output_intermediate_maps+'/meta_{0}/'.format(config['f'])

    # original path for particle counts
    '''
    for the original, let's download the 2048 and downgrade them to 1024.
    shell_ = hp.ud_grade(shell_,nside_out=config['nside_intermediate'])
    '''
    shell = np.load(path_sims+'/run_{0}//shells_nside=512.npz'.format(config['f']))


    from ekit import paths as path_tools
    z_bounds     = dict()
    z_bounds['z-high'] = np.array([shell['shell_info'][i][3] for i in range(len(shell['shell_info']))])
    z_bounds['z-low'] = np.array([shell['shell_info'][i][2] for i in range(len(shell['shell_info']))])

    i_sprt = np.argsort(z_bounds['z-low'])
    z_bounds['z-low']= (z_bounds['z-low'])[i_sprt]
    z_bounds['z-high']= (z_bounds['z-high'])[i_sprt]


    z_bin_edges = np.hstack([z_bounds['z-low'],z_bounds['z-high'][-1]])
    
    '''
    for s_ in (range(len(z_bounds['z-high']))):
        path_ = base+'/lens_{0}_{1}.fits'.format(s_,config['nside_out'])
        if not os.path.exists(path_):
            try:
                shell_ = shell['shells'][i_sprt[s_]]
                shell_ =  (shell_-np.mean(shell_))/np.mean(shell_)
                #shell_ = hp.ud_grade(shell_, nside_out = config['nside_out'])

                fits_f = Table()
                fits_f['T'] = shell_
                if os.path.exists(path_):
                    os.remove(path_)
                fits_f.write(path_)
            except:
                pass
    '''
    path_ = base+'/gg_{0}_{1}.fits'.format(len(z_bounds['z-high'])-1,config['nside_out'])
    if not os.path.exists(path_):
        # SAVE LENS MAPS  *****************************************************************
        for s_ in (range(len(z_bounds['z-high']))):
            path_ = base+'/lens_{0}_{1}.fits'.format(s_,config['nside_out'])
            if not os.path.exists(path_):
                shell_ = shell['shells'][i_sprt[s_]]
                shell_ =  (shell_-np.mean(shell_))/np.mean(shell_)
                shell_ = hp.ud_grade(shell_, nside_out = config['nside_out'])

                fits_f = Table()
                fits_f['T'] = shell_
                if os.path.exists(path_):
                    os.remove(path_)
                fits_f.write(path_)



            path_ = base+'/lens_{0}_{1}.fits'.format(s_,config['nside_intermediate'])
            if not os.path.exists(path_):
                shell_ = shell['shells'][i_sprt[s_]]
                shell_ =  (shell_-np.mean(shell_))/np.mean(shell_)
                shell_ = hp.ud_grade(shell_, nside_out = config['nside_intermediate'])

                fits_f = Table()
                fits_f['T'] = shell_
                if os.path.exists(path_):
                    os.remove(path_)
                fits_f.write(path_)


                
                
                
    # SAVE CONVERGENCE PLANES ********************************************************
    kappa_pref_evaluated = brk.kappa_prefactor(cosmology.H0, cosmology.Om0, length_unit = 'Mpc')
    comoving_edges = [cosmology.comoving_distance(x_) for x_ in np.array((z_bounds['z-low']))]


    z_centre = np.empty((len(comoving_edges)-1))
    for i in range(len(comoving_edges)-1):
        z_centre[i] = z_at_value(cosmology.comoving_distance,0.5*(comoving_edges[i]+comoving_edges[i+1]))

    un_ = comoving_edges[:(i+1)][0].unit
    comoving_edges = np.array([c.value for c in comoving_edges])
    comoving_edges = comoving_edges*un_


    overdensity_array = [np.zeros(hp.nside2npix(config['nside_intermediate']))]


    path_ = base+'/gg_{0}_{1}.fits'.format(len(z_bounds['z-high'])-1,config['nside_out'])
    if not os.path.exists(path_):

        #print ('load lens')
        for s_ in (range(len(z_bounds['z-high']))):
            try:
                path_ = base+'/lens_{0}_{1}.fits'.format(s_,config['nside_intermediate'])
                m_ = pf.open(path_)
                overdensity_array.append(m_[1].data['T'])
            except:
                if shell !=0:
                    overdensity_array.append(np.zeros(hp.nside2npix(config['nside_intermediate'])))
                #pass

        #print ('done ++++')
        overdensity_array = np.array(overdensity_array)
        #print(overdensity_array.shape)

        from bornraytrace import lensing
        kappa_lensing = np.copy(overdensity_array)*0.


        for i in (np.arange(kappa_lensing.shape[0])):
            try:
                kappa_lensing[i] = lensing.raytrace(cosmology.H0, cosmology.Om0,
                                             overdensity_array=overdensity_array[:i].T,
                                             a_centre=1./(1.+z_centre[:i]),
                                             comoving_edges=comoving_edges[:(i+1)])
            except:
                print ('failed kappa ',i)
               

        # make g1 and g2 ---
        for i in (range(kappa_lensing.shape[0])):
            path_ = base+'/gg_{0}_{1}.fits'.format(i,config['nside_out'])
            try:
                 os.remove(path_)
            except:
                pass

            if not os.path.exists(path_):

                g1_, g2_ = gk_inv(kappa_lensing[i]-np.mean(kappa_lensing[i]),kappa_lensing[i]*0.,config['nside_intermediate'],config['nside_intermediate']*2)
                g1_IA, g2_IA = gk_inv(overdensity_array[i]-np.mean(overdensity_array[i]),kappa_lensing[i]*0.,config['nside_intermediate'],config['nside_intermediate']*2)


                fits_f = Table()
                #alm = hp.sphtfunc.map2alm(g1_)
                #g1_ = hp.sphtfunc.alm2map(alm,nside= config['nside_out'])
                fits_f['g1'] = hp.ud_grade(g1_,nside_out =config['nside_out'])
                fits_f['g2'] = hp.ud_grade(g2_,nside_out =config['nside_out'])
                fits_f['g1_IA'] = hp.ud_grade(g1_IA,nside_out =config['nside_out'])
                fits_f['g2_IA'] = hp.ud_grade(g2_IA,nside_out =config['nside_out'])

                fits_f.write(path_)
                path_ = base+'/lens_{0}_{1}.fits'.format(i,config['nside_out'])
                os.system('rm {0}'.format(path_))











    # path where kappa/g1/g2 maps are stored
    base_b = output_intermediate_maps+'/meta_baryons_{0}/'.format(config['f'])

    # original path for particle counts
    shell = np.load(path_sims+'/run_{0}//baryonified_shells.npz'.format(config['f']))


    from ekit import paths as path_tools
    z_bounds     = dict()
    z_bounds['z-high'] = np.array([shell['shell_info'][i][3] for i in range(len(shell['shell_info']))])
    z_bounds['z-low'] = np.array([shell['shell_info'][i][2] for i in range(len(shell['shell_info']))])

    i_sprt = np.argsort(z_bounds['z-low'])
    z_bounds['z-low']= (z_bounds['z-low'])[i_sprt]
    z_bounds['z-high']= (z_bounds['z-high'])[i_sprt]


    z_bin_edges = np.hstack([z_bounds['z-low'],z_bounds['z-high'][-1]])
    path_ = base_b+'/gg_{0}_{1}.fits'.format(len(z_bounds['z-high'])-1,config['nside_out'])

    
    
    if not os.path.exists(path_):
        # SAVE LENS MAPS  *****************************************************************
        for s_ in (range(len(z_bounds['z-high']))):
            path_ = base+'/lens_{0}_{1}.fits'.format(s_,config['nside_out'])
            if not os.path.exists(path_):
                shell_ = shell['shells'][i_sprt[s_]]
                shell_ =  (shell_-np.mean(shell_))/np.mean(shell_)
                shell_ = hp.ud_grade(shell_, nside_out = config['nside_out'])

                fits_f = Table()
                fits_f['T'] = shell_
                if os.path.exists(path_):
                    os.remove(path_)
                fits_f.write(path_)



            path_ = base+'/lens_{0}_{1}.fits'.format(s_,config['nside_intermediate'])
            if not os.path.exists(path_):
                shell_ = shell['shells'][i_sprt[s_]]
                shell_ =  (shell_-np.mean(shell_))/np.mean(shell_)
                shell_ = hp.ud_grade(shell_, nside_out = config['nside_intermediate'])

                fits_f = Table()
                fits_f['T'] = shell_
                if os.path.exists(path_):
                    os.remove(path_)
                fits_f.write(path_)



    # SAVE CONVERGENCE PLANES ********************************************************
    kappa_pref_evaluated = brk.kappa_prefactor(cosmology.H0, cosmology.Om0, length_unit = 'Mpc')
    comoving_edges = [cosmology.comoving_distance(x_) for x_ in np.array((z_bounds['z-low']))]


    z_centre = np.empty((len(comoving_edges)-1))
    for i in range(len(comoving_edges)-1):
        z_centre[i] = z_at_value(cosmology.comoving_distance,0.5*(comoving_edges[i]+comoving_edges[i+1]))

    un_ = comoving_edges[:(i+1)][0].unit
    comoving_edges = np.array([c.value for c in comoving_edges])
    comoving_edges = comoving_edges*un_


    overdensity_array = [np.zeros(hp.nside2npix(config['nside_intermediate']))]


    path_ = base_b+'/gg_{0}_{1}.fits'.format(len(z_bounds['z-high'])-1,config['nside_out'])


    if not os.path.exists(path_):

        #print ('load lens')
        for s_ in (range(len(z_bounds['z-high']))):
            try:
                path_ = base_b+'/lens_{0}_{1}.fits'.format(s_,config['nside_intermediate'])
                m_ = pf.open(path_)
                overdensity_array.append(m_[1].data['T'])
            except:
                if shell !=0:
                    overdensity_array.append(np.zeros(hp.nside2npix(config['nside_intermediate'])))
                #pass






        #print ('done ++++')
        overdensity_array = np.array(overdensity_array)
        #print(overdensity_array.shape)

        from bornraytrace import lensing
        kappa_lensing = np.copy(overdensity_array)*0.


        for i in (np.arange(kappa_lensing.shape[0])):
            try:
                kappa_lensing[i] = lensing.raytrace(cosmology.H0, cosmology.Om0,
                                             overdensity_array=overdensity_array[:i].T,
                                             a_centre=1./(1.+z_centre[:i]),
                                             comoving_edges=comoving_edges[:(i+1)])
            except:
                pass



        # make g1 and g2 ---
        # make g1 and g2 ---
        for i in (range(kappa_lensing.shape[0])):
            path_ = base+'/gg_{0}_{1}.fits'.format(i,config['nside_out'])
            try:
                 os.remove(path_)
            except:
                pass

            if not os.path.exists(path_):

                g1_, g2_ = gk_inv(kappa_lensing[i]-np.mean(kappa_lensing[i]),kappa_lensing[i]*0.,config['nside_intermediate'],config['nside_intermediate']*2)
                g1_IA, g2_IA = gk_inv(overdensity_array[i]-np.mean(overdensity_array[i]),kappa_lensing[i]*0.,config['nside_intermediate'],config['nside_intermediate']*2)


                fits_f = Table()
                #alm = hp.sphtfunc.map2alm(g1_)
                #g1_ = hp.sphtfunc.alm2map(alm,nside= config['nside_out'])
                fits_f['g1'] = hp.ud_grade(g1_,nside_out =config['nside_out'])
                fits_f['g2'] = hp.ud_grade(g2_,nside_out =config['nside_out'])
                fits_f['g1_IA'] = hp.ud_grade(g1_IA,nside_out =config['nside_out'])
                fits_f['g2_IA'] = hp.ud_grade(g2_IA,nside_out =config['nside_out'])

                fits_f.write(path_)
                path_ = base+'/lens_{0}_{1}.fits'.format(i,config['nside_out'])
                os.system('rm {0}'.format(path_))
                
                

    # read n(z) ********************
    mu = pf.open(config['2PT_FILE'])
    #random_rel = np.random.randint(0,6000,1)[0]
    random_rel = 0
    redshift_distributions_sources = {'z':None,'bins':dict()}
    redshift_distributions_sources['z'] = mu[8+random_rel].data['Z_MID']
    for ix in config['sources_bins']:
        redshift_distributions_sources['bins'][ix] = mu[8+random_rel].data['BIN{0}'.format(ix)]
    mu = None

    # prepare g1,g2 maps, interpolate n(z) at the shells location and apply redshift shifts
    g1_tomo = dict()
    g2_tomo = dict()
    d_tomo = dict()
    
    g1_tomo_b = dict()
    g2_tomo_b = dict()
    d_tomo_b = dict()
    
    nz_kernel_sample_dict = dict()
    for tomo_bin in config['sources_bins']:
        g1_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
        g2_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
        d_tomo[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
        
        g1_tomo_b[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
        g2_tomo_b[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
        d_tomo_b[tomo_bin] = np.zeros(hp.nside2npix(config['nside']))
        
        redshift_distributions_sources['bins'][tomo_bin][250:] = 0.
        nz_sample = brk.recentre_nz(np.array(z_bin_edges).astype('float'),  redshift_distributions_sources['z']+config['dz_sources'][tomo_bin-1],  redshift_distributions_sources['bins'][tomo_bin] )
        nz_kernel_sample_dict[tomo_bin] = nz_sample*(z_bin_edges[1:]-z_bin_edges[:-1])

    # these quantities are for the IA computation
    c1 = (5e-14 * (u.Mpc**3.)/(u.solMass * u.littleh**2) ) 
    c1_cgs = (c1* ((u.littleh/(cosmology.H0.value/100))**2.)).cgs
    rho_c1 = (c1_cgs*cosmology.critical_density(0)).value


    # fill in g1,g1 maps
    for i in (range(2,len(comoving_edges)-1)):

            path_ = base+'/lens_{0}_{1}.fits'.format(i,config['nside_out'])
            pathgg_ = base+'/gg_{0}_{1}.fits'.format(i,config['nside_out'])
  
            k_ = pf.open(pathgg_)
            d_ = pf.open(path_)
            IA_f = iaa.F_nla(z_centre[i], cosmology.Om0, rho_c1=rho_c1,A_ia = config['A_IA'], eta=config['eta_IA'], z0=config['z0_IA'],  lbar=0., l0=1e-9, beta=0.)
            #print ((k_[1].data['T']))
            for tomo_bin in config['sources_bins']:         
                m_ = 1.+config['m_sources'][tomo_bin-1]
                b = 1.#extra_params['db'][tomo_bin]
                if SC:
                    g1_tomo[tomo_bin]  +=  ((1.+(b*d_[1].data['T']))*(k_[1].data['g1']+k_[1].data['g1_IA']*IA_f))*nz_kernel_sample_dict[tomo_bin][i]
                    g2_tomo[tomo_bin]  +=  ((1.+(b*d_[1].data['T']))*(k_[1].data['g2']+k_[1].data['g2_IA']*IA_f))*nz_kernel_sample_dict[tomo_bin][i]
                    d_tomo[tomo_bin] +=  (1.+b*d_[1].data['T'])*nz_kernel_sample_dict[tomo_bin][i]
                else:
                    g1_tomo[tomo_bin]  +=  (k_[1].data['g1']+k_[1].data['g1_IA']*IA_f)*nz_kernel_sample_dict[tomo_bin][i]
                    g2_tomo[tomo_bin]  +=  (k_[1].data['g2']+k_[1].data['g2_IA']*IA_f)*nz_kernel_sample_dict[tomo_bin][i]
                    d_tomo[tomo_bin] +=  (1.+d_[1].data['T'])*nz_kernel_sample_dict[tomo_bin][i]


                    
            path_ = base_b+'/lens_{0}_{1}.fits'.format(i,config['nside_out'])
            pathgg_ = base_b+'/gg_{0}_{1}.fits'.format(i,config['nside_out'])
  
            k_ = pf.open(pathgg_)
            d_ = pf.open(path_)
            IA_f = iaa.F_nla(z_centre[i], cosmology.Om0, rho_c1=rho_c1,A_ia = config['A_IA'], eta=config['eta_IA'], z0=config['z0_IA'],  lbar=0., l0=1e-9, beta=0.)
            #print ((k_[1].data['T']))
            for tomo_bin in config['sources_bins']:         
                m_ = 1.+config['m_sources'][tomo_bin-1]
                b = 1.#extra_params['db'][tomo_bin]
                if SC:
                    g1_tomo_b[tomo_bin]  +=  ((1.+(b*d_[1].data['T']))*(k_[1].data['g1']+k_[1].data['g1_IA']*IA_f))*nz_kernel_sample_dict[tomo_bin][i]
                    g2_tomo_b[tomo_bin]  +=  ((1.+(b*d_[1].data['T']))*(k_[1].data['g2']+k_[1].data['g2_IA']*IA_f))*nz_kernel_sample_dict[tomo_bin][i]
                    d_tomo_b[tomo_bin] +=  (1.+b*d_[1].data['T'])*nz_kernel_sample_dict[tomo_bin][i]
                else:
                    g1_tomo_b[tomo_bin]  +=  (k_[1].data['g1']+k_[1].data['g1_IA']*IA_f)*nz_kernel_sample_dict[tomo_bin][i]
                    g2_tomo_b[tomo_bin]  +=  (k_[1].data['g2']+k_[1].data['g2_IA']*IA_f)*nz_kernel_sample_dict[tomo_bin][i]
                    d_tomo_b[tomo_bin] +=  (1.+d_[1].data['T'])*nz_kernel_sample_dict[tomo_bin][i]

                    
                    
                    
    # apply rotations and change nside
    for tomo_bin in config['sources_bins']:
    
        g1_tomo[tomo_bin] = hp.pixelfunc.ud_grade(g1_tomo[tomo_bin],nside_out=nside_out)
        g2_tomo[tomo_bin] = hp.pixelfunc.ud_grade(g2_tomo[tomo_bin],nside_out=nside_out)                   
        d_tomo[tomo_bin] =  hp.pixelfunc.ud_grade(d_tomo[tomo_bin] ,nside_out=nside_out)


   # Add noise ++++++++++++++++++++++++++++++++=
 
    sources_maps = dict()
    for tomo_bin in config['sources_bins']:   
        if noise_type == 'random_depth':
            depth_weigth = np.load('/global/cfs/cdirs/des/mass_maps/Maps_final/depth_maps_Y3_{0}_numbdensity.npy'.format(nside_out),allow_pickle=True).item()
            mcal_catalog = load_obj('/global/cfs/cdirs/des/mass_maps/Maps_final/data_catalogs_weighted_{0}'.format(tomo_bin-1))


            pix_ = convert_to_pix_coord(mcal_catalog['ra'], mcal_catalog['dec'], nside=nside_out)
            
            dp_ = copy.deepcopy(depth_weigth[tomo_bin-1])
            if rot ==1:
                dp_ = rotate_map_approx(depth_weigth[tomo_bin-1],[ 180 ,0 , 0], flip=False,nside =nside_out )
                rot_angles = [180, 0, 0]
                flip=False
                rotU = hp.rotator.Rotator(rot=rot_angles, deg=True)
                alpha, delta = hp.pix2ang(512, np.arange(hp.nside2npix(512)))
                rot_alpha, rot_delta = rotU(alpha, delta)
                if not flip:
                    rot_i = hp.ang2pix(512, rot_alpha, rot_delta)
                else:
                    rot_i = hp.ang2pix(512, np.pi-rot_alpha, rot_delta)
                pix_ = rot_i[pix_]
            if rot ==2:
                dp_ = rotate_map_approx(depth_weigth[tomo_bin-1],[ 90 ,0 , 0], flip=True,nside = nside_out )

                rot_angles = [90, 0, 0]
                flip=True
                rotU = hp.rotator.Rotator(rot=rot_angles, deg=True)
                alpha, delta = hp.pix2ang(512, np.arange(hp.nside2npix(512)))
                rot_alpha, rot_delta = rotU(alpha, delta)
                if not flip:
                    rot_i = hp.ang2pix(512, rot_alpha, rot_delta)
                else:
                    rot_i = hp.ang2pix(512, np.pi-rot_alpha, rot_delta)
                pix_ = rot_i[pix_]
            if rot ==3:
                dp_ = rotate_map_approx(depth_weigth[tomo_bin-1],[ 270 ,0 , 0], flip=True,nside = nside_out)

                rot_angles = [270, 0, 0]
                flip=True
                rotU = hp.rotator.Rotator(rot=rot_angles, deg=True)
                alpha, delta = hp.pix2ang(512, np.arange(hp.nside2npix(512)))
                rot_alpha, rot_delta = rotU(alpha, delta)
                if not flip:
                    rot_i = hp.ang2pix(512, rot_alpha, rot_delta)
                else:
                    rot_i = hp.ang2pix(512, np.pi-rot_alpha, rot_delta)
                pix_ = rot_i[pix_]     
                    
                    
            mask = np.in1d(np.arange(hp.nside2npix(nside_out)),pix_)


     
           # The following bit samples galaxies from the <n> map we just loaded; then , it associates weights equal to the weight of one
       # of the catalogue galaxies in the same pixel . Last, it associates ellipticities based on w-ellipticity relation learned from the catalog

           
            df2 = pd.DataFrame(data = {'w':mcal_catalog['w'] ,'pix_':pix_},index = pix_)
            nn = np.random.poisson(dp_)
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

        
        elif noise_type == 'desy3':
   
            mcal_catalog = load_obj('/global/cfs/cdirs/des/mass_maps/Maps_final/data_catalogs_weighted_{0}'.format(tomo_bin))
            dec1 = mcal_catalog['dec']
            ra1 = mcal_catalog['ra']
            e1 = mcal_catalog['e1']
            e2 = mcal_catalog['e2']
            w = mcal_catalog['w'] 
            pix_ = convert_to_pix_coord(ra1,dec1, nside=nside_out)
            
            if rot ==0:
                pix = copy.deepcopy(pix_)
                
            if rot ==1:
                dp_ = rotate_map_approx(depth_weigth[tomo_bin-1],[ 180 ,0 , 0], flip=False,nside =nside_out )
                rot_angles = [180, 0, 0]
                flip=False
                rot = hp.rotator.Rotator(rot=rot_angles, deg=True)
                alpha, delta = hp.pix2ang(512, pix_)
                rot_alpha, rot_delta = rot(alpha, delta)
                if not flip:
                    pix = hp.ang2pix(512, rot_alpha, rot_delta)
                else:
                    pix = hp.ang2pix(512, np.pi-rot_alpha, rot_delta)
              
         
            if rot ==2:
                dp_ = rotate_map_approx(depth_weigth[tomo_bin-1],[ 90 ,0 , 0], flip=True,nside = nside_out )

                rot_angles = [90, 0, 0]
                flip=True
                rot = hp.rotator.Rotator(rot=rot_angles, deg=True)
                alpha, delta = hp.pix2ang(512, pix_)
                rot_alpha, rot_delta = rot(alpha, delta)
                if not flip:
                    pix = hp.ang2pix(512, rot_alpha, rot_delta)
                else:
                    pix = hp.ang2pix(512, np.pi-rot_alpha, rot_delta)
        
            if rot ==3:
                dp_ = rotate_map_approx(depth_weigth[tomo_bin-1],[ 270 ,0 , 0], flip=True,nside = nside_out)

                rot_angles = [270, 0, 0]
                flip=True
                rot = hp.rotator.Rotator(rot=rot_angles, deg=True)
                alpha, delta = hp.pix2ang(512, pix_)
                rot_alpha, rot_delta = rot(alpha, delta)
                if not flip:
                    pix = hp.ang2pix(512, rot_alpha, rot_delta)
                else:
                    pix = hp.ang2pix(512, np.pi-rot_alpha, rot_delta)
                
            
            del mcal_catalog
            gc.collect() 
            
            
                
            
        unique_pix, idx, idx_rep = np.unique(pix, return_index=True, return_inverse=True)

        # let us sample the noiseless shear maps at the galaxy location
        g1_ = g1_tomo[tomo_bin][pix]
        g2_ = g2_tomo[tomo_bin][pix]
            
        # SC correction factor
        if SC:
            f = 1./np.sqrt(d_tomo[tomo_bin]/np.sum(nz_kernel_sample_dict[tomo_bin]))
            f = f[pix]
        else:
            f = 1.
            
        sources_maps[tomo_bin] = dict()

        n_map_sc = np.zeros(hp.nside2npix(nside_out))
        n_map_sc[unique_pix] += np.bincount(idx_rep, weights=w/f**2)



        es1,es2 = apply_random_rotation(e1/f, e2/f)
        es1a,es2a = apply_random_rotation(e1/f, e2/f)

        x1_sc,x2_sc = addSourceEllipticity({'shear1':g1_,'shear2':g2_},{'e1':es1,'e2':es2},es_colnames=("e1","e2"))

        e1_map = np.zeros(hp.nside2npix(nside_out))
        e2_map = np.zeros(hp.nside2npix(nside_out))
        e1r_map = np.zeros(hp.nside2npix(nside_out))
        e2r_map = np.zeros(hp.nside2npix(nside_out))


        e1_map[unique_pix] += np.bincount(idx_rep, weights= x1_sc*w)
        e2_map[unique_pix] += np.bincount(idx_rep, weights= x2_sc*w)
        e1r_map[unique_pix] += np.bincount(idx_rep, weights=es1a*w)
        e2r_map[unique_pix] += np.bincount(idx_rep, weights=es2a*w)

        mask_sims = n_map_sc != 0.
        e1_map[mask_sims]  = e1_map[mask_sims]/(n_map_sc[mask_sims])
        e2_map[mask_sims] =  e2_map[mask_sims]/(n_map_sc[mask_sims])
        e1r_map[mask_sims]  = e1r_map[mask_sims]/(n_map_sc[mask_sims])
        e2r_map[mask_sims] =  e2r_map[mask_sims]/(n_map_sc[mask_sims])



        m_ = 1.+config['m_sources'][tomo_bin-1]


        print ( m_,tomo_bin)
        sources_maps[tomo_bin]  = {'e1':m_*e1_map,'e2':m_*e2_map,'e1r':m_*e1r_map,'e2r':m_*e2r_map} 

        

        
        
        
        # let us sample the noiseless shear maps at the galaxy location
        g1_ = g1_tomo_b[tomo_bin][pix]
        g2_ = g2_tomo_b[tomo_bin][pix]
            
        # SC correction factor
        if SC:
            f = 1./np.sqrt(d_tomo_b[tomo_bin]/np.sum(nz_kernel_sample_dict[tomo_bin]))
            f = f[pix]
        else:
            f = 1.


        n_map_sc = np.zeros(hp.nside2npix(nside_out))
        n_map_sc[unique_pix] += np.bincount(idx_rep, weights=w/f**2)


        x1_sc,x2_sc = addSourceEllipticity({'shear1':g1_,'shear2':g2_},{'e1':es1,'e2':es2},es_colnames=("e1","e2"))

        e1_map = np.zeros(hp.nside2npix(nside_out))
        e2_map = np.zeros(hp.nside2npix(nside_out))
        e1r_map = np.zeros(hp.nside2npix(nside_out))
        e2r_map = np.zeros(hp.nside2npix(nside_out))


        e1_map[unique_pix] += np.bincount(idx_rep, weights= x1_sc*w)
        e2_map[unique_pix] += np.bincount(idx_rep, weights= x2_sc*w)

        mask_sims = n_map_sc != 0.
        e1_map[mask_sims]  = e1_map[mask_sims]/(n_map_sc[mask_sims])
        e2_map[mask_sims] =  e2_map[mask_sims]/(n_map_sc[mask_sims])



        sources_maps[tomo_bin]['e1b'] = m_*e1_map
        sources_maps[tomo_bin]['e2b'] = m_*e2_map
        

    # save output 
    save_obj(output_temp+p,sources_maps)


   #srun --nodes=4 --tasks-per-node=2  python run_cosmogrid_baryons.py

#salloc --nodes 4 --qos interactive --time 04:00:00 --constraint cpu --account=des



# some config
nside = 512 #nside cosmogrid particle count maps
nside_out = 512 #nside final noisy maps
SC = True #apply SC or not
noise_rels = 10 # number of noise realisations considered
rot_num = 4 # number of rotations considered (max 4)
A_IA = 0.0
e_IA = 0.0
runs_cosmo = 7 # number of cosmogrid independent maps
noise_type = 'desy3' # or 'random_depth'


# this is the path to the cosmogrid sims. Chose among the ones below
#path_sims = '/global/cfs/cdirs/des/cosmogrid/raw/fiducial/cosmo_delta_s8_p/'
#path_sims = '/global/cfs/cdirs/des/cosmogrid/raw/fiducial/cosmo_delta_s8_m/'
#path_sims = '/global/cfs/cdirs/des/cosmogrid/raw/fiducial/cosmo_delta_Om_p/'
#path_sims = '/global/cfs/cdirs/des/cosmogrid/raw/fiducial/cosmo_delta_Om_m/'
#path_sims = '/global/cfs/cdirs/des/cosmogrid/raw/fiducial/cosmo_fiducial/' # this is the fiducal
path_sims = '/global/cfs/cdirs/des/cosmogrid/raw/grid/cosmo_002846/'
# this is the path to intermediate products like kappa g1 g2 maps (already run)
# chose one among the ones below
#output_intermediate_maps = '/global/cfs/cdirs/des/mgatti/cosmogrid/cosmo_delta_s8_m/'
#output_intermediate_maps = '/global/cfs/cdirs/des/mgatti/cosmogrid/cosmo_delta_s8_p/'
#output_intermediate_maps = '/global/cfs/cdirs/des/mgatti/cosmogrid/cosmo_delta_Om_p/'
#output_intermediate_maps = '/global/cfs/cdirs/des/mgatti/cosmogrid/cosmo_delta_Om_m/'
output_intermediate_maps = '/global/cfs/cdirs/des/mgatti/cosmogrid_002846_1024_3/' # this is the fiducial run


# the final noisy maps will be saved here

output_temp = '/global/cfs/cdirs/des/mgatti/cosmogrid/baryons_002846_final/'

if not os.path.exists(output_intermediate_maps):
    try:
        os.mkdir(output_intermediate_maps)
    except:
        pass




if __name__ == '__main__':



    import glob
    runstodo=[]
    count = 0
    miss = 0



    for f in range(runs_cosmo):

        if not os.path.exists(output_intermediate_maps+'/meta_{0}/'.format(f)):
            try:
                os.mkdir(output_intermediate_maps+'/meta_{0}/'.format(f))
            except:
                pass

        if not os.path.exists(output_intermediate_maps+'/meta_baryons_{0}/'.format(f)):
            try:
                os.mkdir(output_intermediate_maps+'/meta_baryons_{0}/'.format(f))
            except:
                pass


        with open(path_sims+'/run_{0}/params.yml'.format(f), "r") as f_in:
            config = yaml.safe_load(f_in.read())

        Omegam = config['Om']
        s8 = config['s8']
        ns = config['ns']
        Ob = config['Ob']
        h = config['H0']
        w0 = config['w0']



        for i in range(rot_num):

            for nn in range(noise_rels):

                params_dict = dict()
                params_dict['Omegam'] = np.float(Omegam)
                params_dict['s8'] = np.float(s8)

                params_dict['A'] = A_IA
                params_dict['E'] = e_IA
                params_dict['noise'] = nn


                params_dict['rot'] = i
                params_dict['ns'] = np.float(ns)
                params_dict['h'] = np.float(h)
                params_dict['ob'] = np.float(Ob)

                params_dict['SC'] = SC
                params_dict['f'] = f


                params_dict['m1'] =  -0.002
                params_dict['m2'] = -0.017
                params_dict['m3'] = -0.029
                params_dict['m4'] = -0.038

                params_dict['dz1'] = 0.
                params_dict['dz2'] = 0.
                params_dict['dz3'] = 0.
                params_dict['dz4'] = 0.


                p = str(f)+'_'+noise_type+'_'+str(Omegam )+'_'+str(s8)+'_'+str(ns)+'_'+str(Ob)+'_'+str(h )+'_'+str(A_IA )+'_'+str(e_IA )+'_w'+str(w0)+'_'+str(i+1)+'_noise_'+str(nn)


                if not os.path.exists(output_temp+p+'.pkl'):
                    runstodo.append([p,params_dict])
                    miss+=1
                else:
                    count +=1





    print (len(runstodo),count,miss)


    run_count=0
    
    from mpi4py import MPI
    while run_count<len(runstodo):
        comm = MPI.COMM_WORLD
#
        if (run_count+comm.rank)<len(runstodo):
            #try:
                make_maps(runstodo[run_count+comm.rank])
           # except:
           #     pass
        #if (run_count)<len(runstodo):
        #    make_maps(runstodo[run_count])
        run_count+=comm.size
        comm.bcast(run_count,root = 0)
        comm.Barrier()
       

   # while run_count<len(runstodo):
   #     
##
   #     if (run_count)<len(runstodo):
   #         #try:
   #             make_maps(runstodo[run_count])
   #        # except:
   #        #     pass
   #     #if (run_count)<len(runstodo):
   #     #    make_maps(runstodo[run_count])
   #     run_count+=1
   #     #comm.bcast(run_count,root = 0)
   #     #comm.Barrier()
##srun --nodes=4 --tasks-per-node=7  python run_cosmogrid_baryons.py
##srun --nodes=1 --tasks-per-node=4 --cpus-per-task=16 --cpu-bind=cores  python run_cosmogrid_baryons.py

