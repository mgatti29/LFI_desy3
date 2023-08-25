def write_values_sp(path, params, idd,neutrino):
    a = open(path+'/values_{0}.ini'.format(idd),'w')
    a.write('[cosmological_parameters]\n')
    a.write('omega_m = {0}\n'.format(params['om']))
    a.write('sigma8_input = {0}\n'.format(params['s8']))
    a.write('h0 = {0}\n'.format(params['h']))
    a.write('omega_b = {0}\n'.format(params['ob']))
    a.write('n_s = {0}\n'.format(params['ns']))
    a.write('w = {0}\n'.format(params['w']))
    a.write('A_s = {0}\n'.format(2.19e-9))
    #a.write('omnuh2 = {0}\n'.format(neutrino))
    #a.write('massive_nu = {0}\n'.format(3))
    
    
    
    a.write('mnu = {0}\n'.format(neutrino))
    a.write('massive_nu = {0}\n'.format(3))
    a.write('nnu = {0}\n'.format(3.046))
    
    ## mnu = 0.02794    0.056    0.0745
    ## num_massive_neutrinos = 1
    ## nnu = 1.046



    a.write('omega_k = {0}\n'.format(0.))
    
    
    a.write('tau = {0}\n'.format(0.0697186))
    a.write('yhe = {0}\n'.format(0.245341))
         
    a.write('[shear_calibration_parameters]\n')
    for j in range(4):
        a.write('m{0} = {1}\n'.format(j+1,1))
           
    a.write('[wl_photoz_errors]\n')
    for j in range(4):
        a.write('bias_{0} = {1}\n'.format(j+1,0.))
                     
            
    a.write('[bias_lens]\n')
    for j in range(6):
        a.write('b{0} = {1}\n'.format(j+1,np.random.uniform( 0.  ,3)))
       

        
    a.write('[mag_alpha_lens]\n')
    a.write('alpha_1 = 1.21433242 \n')
    a.write('alpha_2 = 1.1485733\n')
    a.write('alpha_3 = 1.87589872\n')
    a.write('alpha_4 = 1.96942553\n') 
    a.write('alpha_5 = 1.78050552\n') 
    a.write('alpha_6 = 2.47892785\n')

    a.write('[intrinsic_alignment_parameters]\n')
    a.write('z_piv = 0.62 \n')

    
    a.write('A1 = {0}\n'.format(0.))
    a.write('A2 = {0}\n'.format(0.))
    a.write('alpha1 = {0}\n'.format(0.))
    a.write('alpha2 = {0}\n'.format(0.))
    a.write('bias_ta = {0}\n'.format(0.))
    a.write('bias_tt = {0}\n'.format(0.))
    a.write('Adel = {0}\n'.format(0.))
    
   
   

    a.write('[lens_photoz_errors]\n')
    a.write('bias_1 =  {0}\n'.format(np.random.normal(-0.009, 0.007)))
    a.write('bias_2 =  {0}\n'.format(np.random.normal(-0.035, 0.011)))
    a.write('bias_3 =  {0}\n'.format(np.random.normal(-0.005, 0.006)))
    a.write('bias_4 =  {0}\n'.format(np.random.normal(-0.007, 0.006)))
    a.write('bias_5 =  {0}\n'.format(np.random.normal(0.002 ,0.007 )))
    a.write('bias_6 =  {0}\n'.format(np.random.normal(0.002 ,0.008 )))
    a.write('width_1 = {0}\n'.format(np.random.normal(0.975 , 0.062)))
    a.write('width_2 = {0}\n'.format(np.random.normal(1.306 , 0.093)))
    a.write('width_3 = {0}\n'.format(np.random.normal(0.87  ,0.054 )))
    a.write('width_4 = {0}\n'.format(np.random.normal(0.918 , 0.051)))
    a.write('width_5 = {0}\n'.format(np.random.normal(1.08  ,0.067 )))
    a.write('width_6 = {0}\n'.format(np.random.normal(0.845 , 0.073)))
    
    a.close()
    
            
#def write_params_sp(path,  idd,rel,folder = 'test'):
#    a = open(path+'/params_SR_{0}.ini'.format(idd),'w')
#    a.write('%include /global/u2/m/mgatti/Mass_Mapping/Moments_analysis/mcmc_PWH/Y3_params_priors/params.ini\n')
#    a.write('[fits_nz]\n')
#    a.write('file = cosmosis-standard-library/number_density/load_nz_fits/load_nz_fits_hyp.py\n')
#    a.write('nz_file = %(2PT_FILE)s\n')
#    a.write('data_sets = lens source\n')
#    a.write("rel = '{0}'\n".format(rel))
#    a.write('prefix_section = T\n')
#    a.write('prefix_extension = T\n')
#    a.write('[test]\n')
#    a.write('save_dir={0}\n'.format(folder))
#
#
#
#    a.close()
# 


            
            
            
            
     
    
def write_params_sp(path,  idd,rel,folder = 'test'):
    a = open(path+'/params_SR_{0}.ini'.format(idd),'w')
    a.write('%include /global/cfs/cdirs/des/mxlin/temp/Y3_params_priors/params.ini\n')
    a.write('[fits_nz]\n')
    a.write('file =  /global/cfs/cdirs/des/mxlin/temp/cosmosis-standard-library/number_density/load_nz_fits/load_nz_fits.py\n') #load_nz_fits_hyp
    a.write('nz_file = %(2PT_FILE)s\n')
    a.write('data_sets = lens source\n')
    a.write("rel = '{0}'\n".format(rel))
    a.write('prefix_section = T\n')
    a.write('prefix_extension = T\n')
    a.write('[test]\n')
    a.write('save_dir={0}\n'.format(folder))



    a.close()
 


            
            
            
            
            
            
            

#srun --nodes=1 --tasks-per-node=64 --cpus-per-task=1 --cpu-bind=cores   python test_output_LFI_final.py

#srun --nodes=1 --tasks-per-node=64 --cpus-per-task=1 --cpu-bind=cores   python test_output_LFI_cosmogrid.py

import pickle
def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        f.close()
        
def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        mute =  pickle.load(f)
        f.close()
    return mute



    
   

    
    
    
import numpy as np
import pickle
from pyDOE import *
import os
import sys
sys.path.append('./auxiliary_routines/')
#from routines import *

import glob
'''
This code runs cosmosis to evaluate the DV at different training points.
1 ) First check auxiliary_routines/routines/; the function ``write_values'' takes the parameters value from one of the latin hypercube entry and fills up a values.ini file for cosmosis with those parameters, making sure missing parameters are filled with standard value. It should work out of the box but you might want to check all the parameters you care are there. 
2) 

'''
    
# this folder will be used to store the cosmosis predictions.    
folder_output = '/global/cfs/cdirs/des/mgatti/Dirac/output_SR_new//test'  
    
# path to the emu_cosmosis folder.
path = '/global/cfs/cdirs/des/mxlin/temp/Y3_params_priors/'
output_folder_save = '/global/cfs/cdirs/des/mgatti/Dirac/output_folder_cl_predictions_Dirac_sept23_2_noz25/'
if __name__ == '__main__':
    
    # set your own. Only datafile is really important, priors and scales could be whatever.
    os.environ['DATAFILE']='/global/cfs/cdirs/des//www/y3_chains/data_vectors/2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate_6000HR_mod.fits'
    os.environ['DATAFILE_SR']='/global/u2/m/mgatti/Mass_Mapping/Moments_analysis/mcmc_PWH/data/2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate_sr.npy'
    os.environ['INCLUDEFILE'] = path+'params.ini'
    os.environ['PRIORSINCLUDE']=path+'priors.ini'


    files_ = glob.glob('/global/cfs/cdirs/des/mgatti/Dirac/output_moments_new_dirac_C///*')
    files = []
    for f in files_:
        if 'noiserel30' in f and 'rel0.' in f:
            files.append(f)
            
    run_count = 0
    f = 0
    training_points = len(files)
         
    print (training_points)
    
    
    
    

    run_count = 0


    
    #'''
    from mpi4py import MPI 
    while run_count<training_points:
   
        comm = MPI.COMM_WORLD
        if run_count+comm.rank<(training_points):
        #'''
        #if 1==1:

            #idx =  375#run_count+comm.rank
            idx = run_count+comm.rank


            #m = load_obj(files[idx].split('.pkl')[0])
            label = files[idx].split('/')[-1].split('.pkl')[0]
            fold = label.split('runs')[1].split('_')[0]
            run_= int(label.split('run')[2].split('_')[0])

            print (fold,run_)

            f = open(('/global/homes/m/mgatti/Mass_Mapping/peaks/params_run_1_Niall_{0}.txt'.format(fold)),'r')
            mv_ = []
            h_ = []
            ob_ = []
            ns_ = []
            w_ = []
            om_ =[]
            sigma8_ = []
            for i,f_ in enumerate(f):
                if i>0:
                    ns_.append(float(f_.split(',')[5]))

                    h_.append(float(f_.split(',')[4]))
                    om_.append(float(f_.split(',')[0]))
                    ob_.append(float(f_.split(',')[3]))
                    w_.append(float(f_.split(',')[2]))
                    sigma8_.append(float(f_.split(',')[1]))
                    try:
                        mv_.append(float(f_.split(',')[6]))
                    except:
                        mv_.append(0.06)




            mv = mv_[run_-1]
            h = h_[run_-1]
            neutrino = mv#*h*h/94  

            params = dict()
            params['om'] = om_[run_-1]
            params['s8'] = sigma8_[run_-1]
            params['h'] = h_[run_-1]
            params['ob'] = ob_[run_-1]
            params['ns'] = ns_[run_-1]
            params['w'] = w_[run_-1]



            write_values_sp(path,params,idx,neutrino)

            write_params_sp(path, idx,0,output_folder_save+label+'/')

            os.system('cosmosis Y3_params_priors/params_SR_{1}.ini -p  pipeline.values="{2}/values_{1}.ini"'.format(folder_output+str(label)+'1n',idx,path))


            os.remove('{0}/values_{1}.ini'.format(path,idx))
            os.remove('{0}/params_SR_{1}.ini'.format(path,idx))



        run_count+=comm.size
        comm.bcast(run_count,root = 0)
        comm.Barrier() 
        
    '''
    os.remove('{0}/values_{1}.ini'.format(path,idx))
    os.remove('{0}/params_SR_{1}.ini'.format(path,idx))
    
    
    #os.system('cosmosis Y3_params_priors/params_SR_{1}.ini -p  pipeline.values="{2}/values_{1}.ini"'.format(folder_output+str(label),iii,path))
    '''
    '''
        
 


        
    fail = []
    from mpi4py import MPI 
    while run_count<training_points:
    #while run_count<5: #len(hypercube_samples):
        comm = MPI.COMM_WORLD
        if run_count+comm.rank<(training_points):
        #if run_count<len(hypercube_samples):
            #try:
                m = load_obj(files[run_count+comm.rank].split('.pkl')[0])

                label = files[run_count+comm.rank].split('/')[-1].split('.pkl')[0]

                fold = label.split('runs')[1].split('_')[0]
                
                if (fold !='C') and (fold !='E'):
                    run_= np.int(label.split('run')[2].split('_')[0])

                    passe = False
                    try:
                        f = open(('/global/homes/m/mgatti/Mass_Mapping/peaks/params_run_1_Niall_{0}.txt'.format(fold)),'r')
                        mv_ = []
                        h_ = []
                        for i,f_ in enumerate(f):
                            if i>0:
                                mv_.append(float(f_.split(',')[6]))
                                h_.append(float(f_.split(',')[4]))
                        mv = mv_[run_]
                        h = h_[run_]
                        neutrino = mv*h*h/94  
                        passe = True
                    except:
                        print ('failed folder ',fold)

                    if passe:
                        write_values_sp(path,m[1],run_count+comm.rank,neutrino)
                        write_params_sp(path, run_count+comm.rank,m[1]['hyperrank_rel'],output_folder_save+label+'/')
                        os.system('cosmosis Y3_params_priors/params_SR_{1}.ini -p smallratio_like.save_ratio_filename="{0}" pipeline.values="{2}/values_{1}.ini"'.format(folder_output+str(label),run_count+comm.rank,path))
                        os.remove('{0}/values_{1}.ini'.format(path,run_count+comm.rank))
                        os.remove('{0}/params_SR_{1}.ini'.format(path,run_count+comm.rank))
            #except:
            #    pass

        run_count+=comm.size
        comm.bcast(run_count,root = 0)
        comm.Barrier() 

    '''

    
    '''
  
    while run_count<training_points:
        if run_count<training_points:
                
            m = load_obj(files[run_count].split('.pkl')[0])
            
            label = files[run_count].split('/')[-1].split('.pkl')[0]
            #if not os.path.exists(folder_output+str(label)+'.npy'):
            if 1==1:
                write_values_sp(path,m[1],run_count)
                write_params_sp(path, run_count,m[1]['hyperrank_rel'],output_folder_save+label+'/')
                os.system('cosmosis Y3_params_priors/params_SR_{1}.ini -p smallratio_like.save_ratio_filename="{0}" pipeline.values="{2}/values_{1}.ini"'.format(folder_output+str(label),run_count,path))
                os.remove('{0}/values_{1}.ini'.format(path,run_count))
                os.remove('{0}/params_SR_{1}.ini'.format(path,run_count))
        run_count+=1
            
    '''