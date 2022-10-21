#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 13:46:18 2021
@author: tazio
This will calculate the parameters for the Dressler-Shectman test (D-S, 1988)
It can also be used to alculate the deviations of a cluster with randomised velocities
"""

from astropy.table import Table
import numpy as np
import sys, os.path
from decorators import time_it
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.constants import c
c = c.value/1000

def calc_deviation(raT, decT, tbl, n, v_cl, disp_cl ):
    """
    This will make a list with the n nearest galaxies to given coordinates, inclusive
    Input:
        raT (float) - RA of the target
        decT (float) - DEC of the target
            We want to find the closest of all galaxies to these coordinates
        nearest (table) - Table with all the galaxies
        n (int) - Will return the main galaxy and its n nearest neighbors
        v_cl (float) - central velocity of the cluster
                       the derived central velocity will be compared with this
        disp_cl (float) - velocity dispersion of the cluster
                          the derived velocity disp. will be compared with this
    Return:
        d (float) - the deviation of the first galaxy in the list
    """

    new_nearest = tbl.copy()
    curP = SkyCoord(raT, decT, unit='deg')
    new_nearest['distance'] = curP.separation(SkyCoord(tbl['RAn'], tbl['DECn'], unit='deg'))

    # Take the nearest n
    new_nearest.sort(['distance'])
    new_nearest = new_nearest[:n+1]

    # Claculate the velocity dispersion of the group
    disp_gr = do_elmo(new_nearest)
    
    # Calculate the deviation
    d = np.sqrt(((n+1)/(disp_cl**2)) * ((new_nearest['Veln'].mean()-v_cl)**2 + (disp_gr-disp_cl)**2 ))
    return d

def run_original(dc, v_cl, disp_cl, n):
    for i in range(len(dc)):
        dc['deviation'][i] = calc_deviation(dc[i]['RAn'], dc[i]['DECn'], dc, n, v_cl, disp_cl )
    return dc

def run_random( dc, v_cl, disp_cl, n , Nr, start_num ):
    with open(path+'/ds_test/monte_carlo/random_realisations.dat','a') as rr:
        print('Starting with', Nr,'random runs!')
        for num in range(Nr):
            drc = rand_vels( dc,v_cl,disp_cl ) #randomise the velocities
            for i in range(len(drc)):
                drc['deviation'][i] = calc_deviation(drc[i]['RAn'], drc[i]['DECn'], drc, n, v_cl, disp_cl )
    
            # This can be used to write out the randomised table.
            #drc.write(path+'ds_test/monte_carlo_1/run_'+str(num+start_num)+'.fits')
            
            rr.write(str(num+start_num)+'\t'+"%.5f" %(np.sum(drc['deviation'])/len(drc))+'\n')


def get_start_num(out_file): #find the number of the last run
    of = Table.read(out_file,format='ascii.no_header')
    if len(of) == 0:
        start_num=0
    else:
        start_num = (of['col1'][-1])+1
    del of
    return  start_num

def rand_vels(drc, v_cl, disp_cl): 
    """
    Randomise the velocities of the galaxies. 
    For a check of the significance of the substructure
    It will draw the random values from a Normal distribution
    with the parameters of the cluster
    Input:
        drc (table)
        v_cl (float) - central velcity of the cluster
        disp_cl (float) - velcity dispertion of the cluster
    Return:
        drc (table)
    """
    drc['Veln'] = np.random.normal(v_cl,disp_cl,len(dc))
    return(drc)


def do_elmo(dat):
    # Calculate the velocity dispersion as in Tempel et al. 2014
    v_med = np.mean(dat['Veln'])
    dispE = np.sqrt( (np.sum((dat['Veln']-v_med)**2)) / (((1+np.mean(dat[z]))**2)*(len(dat)-1)) )
    return dispE

def check_input(tbl):
    # Run a basic sanity check on the input parameters
    check_ra = (abs(tbl[ra].mean() - Ira_cl) / Ira_cl > 0.1)
    check_dec = (abs(tbl[dec].mean() - Idec_cl) / Idec_cl > 0.1)
    check_vel = (abs(tbl[vel].mean() - Iv_cl) / Iv_cl > 0.1)
    return check_ra, check_dec, check_vel

@time_it
def run_it(dc, ra, dec, vel, Iv_cl, Idisp_cl, Ira_cl, Idec_cl, n, data=True, random=False):
    dc['distance']=0.0
    dc['deviation']=0.0
    # Normalize the cluster parameters
    vel_norm_const = dc[vel].max() - dc[vel].min()
    
    disp_cl = Idisp_cl/vel_norm_const
    v_cl = (Iv_cl - dc[vel].min()) / vel_norm_const
    ra_cl = (Ira_cl - dc[ra].min())/ (dc[ra].max() - dc[ra].min())
    dec_cl = (Idec_cl - dc[dec].min()) / (dc[dec].max() - dc[dec].min())
    
    # Normalise the three axes
    a = dc[ra] - dc[ra].min()
    dc['RAn'] = a/a.max()
    a = dc[dec] - dc[dec].min()
    dc['DECn'] = a/a.max()
    a = dc[vel] - dc[vel].min()
    dc['Veln'] = a/a.max()
    del a

    print('!!!!!!!!!!!!!!!!!!!\nWill be looking for groups with {} members'.format(n))
    
    if data:
        dc = run_original(dc, v_cl, disp_cl, n) 
        dc.write(path+'ds_test/ds_test.mpa_jhu_A1367.N'+str(int(n))+'.fits', overwrite='True')
    if random:
        # Use this to generate a set with randomised velocities and run 
        # the DS test over those
        
        #This is the number of random sets
        Nr = 1#00
        
        #FOR THE MONTE-CARLO: It will append everything at the end of this file, 
        # and the run table name will start with the last number+1
        out_file = path+'ds_test/monte_carlo/random_realisations.dat'
        if not os.path.isfile(out_file):
            with open(path+'ds_test/monte_carlo/random_realisations.dat','w') as rr:
                pass
            start_num = 0
        else:
            start_num = get_start_num(out_file)
        run_random( dc, v_cl, disp_cl, n , Nr, start_num ) 

    return dc


if __name__ == '__main__':
    ################################################
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #   ADJUST THOSE PARAMETERS AS NEEDED
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # Path and table with all sources
    path = '/home/tazio/works/2022/a1367_cluster/'
    dc = Table.read(path+'mpa_jhu_AGES.fits')
    print('{} sources read in'.format(len(dc)))
    
    # The column names containing the needed coordinates
    ra = 'RA'
    dec = 'DEC'
    vel = 'Vel'
    z = 'Z'

    # PARAMETERS OF THE CLUSTER
    # Central velocity km/s from NED
    # this is wrt CMB as it corresponds best to the SDSS velocities
    # ...
    Iv_cl = 6793
    # Velocity dispersion of the cluster km/s from Girardi+'98
    Idisp_cl = 798 
    # Coordinates of the center of the cluster
    Ira_cl = 176.152083
    Idec_cl = 19.758889
    
    # Select only the velocity range of the cluster 
    dc = dc[(dc[vel]>4000) & (dc[vel]<9000)]
    
    # Do you want to run the DS test on the data or
    # you want to run monte carlo simulations
    # Those take the positions on the sky of the real galaxies
    # but assign a random velocities to them, drawn from a Gaussian with
    # the v and disp of the cluster. It then stores the average deviation
    # in a file. Do this many times and you can distinguish the excess deviation
    # !!! This runs ~17 second per one cluster, beit real or randomised.
    data = True
    random = False
    
    # Number of galaxies in each group. Adjust this if needed 
    # though sqrt(N) is a sound value
    n = int(np.sqrt(len(dc))) 

    ################################################
    if any(check_input(dc)):
        proceed = 'n'
        print('\n','!'*20)
        print('Input cluster parameters don\'t seem to match the input table')
        print('(deviation of more than 10% w.r.t. the mean)')
        proceed = input('Should I proceed? ... y/[n] ')
        if proceed != 'y':
            sys.exit()
        else:
            pass

    print('N cluster members =', len(dc))
    print('Working with v_cl =', Iv_cl, 'and disp_cl =', Idisp_cl)
    dc = run_it(dc, ra, dec, vel, Iv_cl, Idisp_cl, Ira_cl, Idec_cl, n, data=data, random=random)

