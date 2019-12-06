#!/wsu/el7/pre-compiled/python/3.7/bin/python

import numpy as np
import copy
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource
from matplotlib import cm
from matplotlib import animation


#plt.switch_backend('agg')
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['xtick.major.width'] = 2.0
mpl.rcParams['ytick.major.width'] = 2.0
mpl.rcParams['xtick.labelsize'] = 13.0
mpl.rcParams['ytick.labelsize'] = 13.0

deg2rad = np.pi / 180.0


class Results() :

    def __init__(self) :

        self.ndir       = 62 #: define me
        self.ntimes     = 0  #: define me 
        self.nocc_alpha = 0  #: define me 
        self.nocc_beta  = 0  #: define me 
        self.loss2   = {}    #: loss from occupied
        self.loss    = {}    #: instantaneous loss
        self.rate    = {}    #: instantaneous rates
        self.results = {}    #: RESULTS*.OUT file
        self.field   = {}    #: field info read in from RESULTS*.OUT
        self.holes   = {}
        self.parts   = {}
        self.summary = {}    #: ground and excited state population
        self.pop = {}        #: population analysis



    def read_loss2( self, rdir='none', edir='none', start=1, end=0 ) :

        """ Read LOSS2-e#-d#.out files """

        if rdir == 'none' :
            print(' ERROR:  need to specify rdir="directory_of_files"' )
            return
        if edir == 'none' :
            print( ' ERROR:  need to specify edir="int" for FILE-e"edir"-d#.out ' )
            return

        if end == 0 : end = self.ndir + 1 

        for idir in range( start, end+1 ) :
            
            myfile = "{fdir}/LOSS2-e{edir}-d{idir}.out".format( fdir=rdir, edir=str(edir), idir=idir )
            iloss2 = {}

            with open( myfile, 'r' ) as infile:
                
                key_heading = infile.readline() 
                key    = [ ikey for ikey in key_heading.split()[:] ]
                values = [ [] for ivalue in range( len(key) ) ]
                iloss2 = dict( zip(key,values) )

                all_lines = infile.readlines()
                for iline in all_lines[:-2] :
                    ilist = iline.split()
                    for i in range(len(key) ) : iloss2[key[i]].append( float(ilist[i]) )
                        
            self.loss2[idir] = copy.deepcopy( iloss2 )
            
            
            
    def read_rates2( self, rdir='none', edir='none', start=1, end=0 ) :

        """ read RATE2-e#-d#.out files """
        
        if rdir == 'none' :
            print(' ERROR:  need to specify rdir="directory_of_files"' )
            return
        if edir == 'none' :
            print( ' ERROR:  need to specify edir="int" for FILE-e"edir"-d#.out ' )
            return

        if end == 0 : end = self.ndir  
            
        for idir in range( start, end+1 ) :

            myfile = "{fdir}/RATE2-e{edir}-d{idir}.out".format( fdir=rdir, edir=str(edir), idir=idir )
            irate = {}

            with open( myfile, 'r' ) as infile:
                
                key_heading = infile.readline() 
                key    = [ ikey for ikey in key_heading.split()[:] ]
                values = [ [] for ivalue in range( len(key) ) ]
                irate  = dict( zip(key,values) )

                all_lines = infile.readlines()
                for iline in all_lines[:] :
                    ilist = iline.split()
                    for i in range(len(key) ) : irate[key[i]].append( float(ilist[i]) )
                        
            self.rate[idir] = copy.deepcopy( irate )



    def read_loss( self, rdir='none', edir='none', start=1, end=0 ) :

        """ read LOSS-e#-d#.out files """

        if rdir == 'none' :
            print(' ERROR:  need to specify rdir="directory_of_files"' )
            return
        if edir == 'none' :
            print( ' ERROR:  need to specify edir="int" for FILE-e"edir"-d#.out ' )
            return
            
        if end == 0 : end = self.ndir  

        for idir in range( start, end+1 ) :

            myfile = "{fdir}/LOSS_PY-e{edir}-d{idir}.out".format( fdir=rdir, edir=str(edir), idir=idir )
            iloss  = {}

            with open( myfile, 'r' ) as infile:
                
                key_heading = infile.readline() 
                key    = [ ikey for ikey in key_heading.split()[:] ]
                values = [ [] for ivalue in range( len(key) ) ]
                iloss = dict( zip(key,values) )

                all_lines = infile.readlines()
                for iline in all_lines :
                    ilist = iline.split()
                    for i in range(len(key) ) : iloss[key[i]].append( float(ilist[i]) )
                        
            self.loss[idir] = copy.deepcopy( iloss )



    def read_summary_ci( self, rdir='none', edir='none', start=1, end=0 ) :

        if rdir == 'none' :
            print(' ERROR:  need to specify rdir="directory_of_files"' )
            return
        if edir == 'none' :
            print( ' ERROR:  need to specify edir="int" for FILE-e"edir"-d#.out ' )
            return

        if end == 0 : end = self.ndir  
            
        for idir in range( start, end+1 ) :

            myfile = "{fdir}/SUMMARY_CI-e{edir}-d{idir}.out".format( fdir=rdir, edir=str(edir), idir=idir )
            isummary = {}

            with open( myfile, 'r' ) as infile :

                key_heading = infile.readline()
                key = [ ikey for ikey in key_heading.split()[:] ]
                values = [ [] for ival in range( len(key) ) ]
                isummary = dict( zip( key, values) )
    
                all_lines = infile.readlines()
                for iline in all_lines :
                    ilist = iline.split()
                    for i in range(len(key)) : isummary[key[i]].append( float(ilist[i]) ) 
                
            self.summary[idir] = copy.deepcopy( isummary )
            


    def read_results( self, rdir='none', edir='none', start=1, end=0 ) :

        if rdir == 'none' :
            print(' ERROR:  need to specify rdir="directory_of_files"' )
            return
        if edir == 'none' :
            print( ' ERROR:  need to specify edir="int" for FILE-e"edir"-d#.out ' )
            return        

        if end == 0 : end = self.ndir  
            
        for idir in range( start, end+1 ) :
               
            myfile = "{fdir}/RESULTS-e{edir}-d{idir}.out".format( fdir=rdir, edir=str(edir), idir=idir )
            iresults = {}

            with open( myfile, 'r' ) as infile :
                
                #: read in field data
                field_line = infile.readline()
                keys0      = [ ikey0 for ikey0 in field_line.split()[1:] ]
                val_line0  = infile.readline()
                values0    = [ float(ival0) for ival0 in val_line0.split()[1:] ]
                ifield = dict( zip( keys0, values0 ) )

                #: read in Results data
                key_line = infile.readline()
                keys   = [ ikey for ikey in key_line.split()[1:] ]
                values = [ [] for i in range( len(keys) ) ]
                iresults = dict( zip(keys,values) )
                
                all_lines = infile.readlines()
                for iline in all_lines :
                    ilist = iline.split()[1:]
                    for i in range(len(values)) : iresults[keys[i]].append( float(ilist[i]) )

            self.results[idir] = copy.deepcopy( iresults )
            self.field[idir] = copy.deepcopy( ifield )


                    
    def read_holes( self, rdir='none', edir='none' , start=1, end=0 ) :

        if rdir == 'none' :
            print(' ERROR:  need to specify rdir="directory_of_files"' )
            return
        if edir == 'none' :
            print( ' ERROR:  need to specify edir="int" for FILE-e"edir"-d#.out ' )
            return

        if end == 0 : end = self.ndir  

        for idir in range( start, end+1 ) :
            
            myfile = "{fdir}/HOLE_POP-e{edir}-d{idir}.out".format( fdir=rdir, edir=str(edir), idir=idir )
            iholes = {}

            with open( myfile, 'r' ) as infile:

                key_heading = infile.readline()
                keys   = key_heading.split()
                values = [ [] for ivalue in range(len(keys)) ]
                iholes = dict( zip(keys,values) )
                
                all_lines = infile.readlines()
                for iline in all_lines :
                    ilist = iline.split()
                    for i in range(len(keys)) : iholes[keys[i]].append( float(ilist[i]) )
            
            self.holes[idir] = copy.deepcopy( iholes )



    def read_parts( self, rdir='none', edir='none', start=1, end=0 ) :

        if rdir == 'none' :
            print(' ERROR:  need to specify rdir="directory_of_files"' )
            return
        if edir == 'none' :
            print( ' ERROR:  need to specify edir="int" for FILE-e"edir"-d#.out ' )
            return

        if end == 0 : end = self.ndir  

        for idir in range( start, end+1 ) :
                        
            myfile = "{fdir}/PART_POP-e{edir}-d{idir}.out".format( fdir=rdir, edir=str(edir), idir=idir )
            iparts = {}
            
            with open( myfile, 'r' ) as infile:

                key_heading = infile.readline()
                keys   = key_heading.split()
                values = [ [] for ivalue in range(len(keys)) ]
                iparts = dict( zip(keys,values) )
                
                all_lines = infile.readlines()
                for iline in all_lines :
                    ilist = iline.split()
                    for i in range(len(keys)) : iparts[keys[i]].append( float(ilist[i]) )
            
            self.parts[idir] = copy.deepcopy( iparts )



    def read_pop( self, rdir='none', edir='none', start=1, end=0 ) :

        if rdir == 'none' :
            print(' ERROR:  need to specify rdir="directory_of_files"' )
            return
        if edir == 'none' :
            print( ' ERROR:  need to specify edir="int" for FILE-e"edir"-d#.out ' )
            return

        if end == 0 : end = self.ndir  

        for idir in range( start, end+1 ) :
                        
            myfile = "{fdir}/POP-e{edir}-d{idir}.out".format( fdir=rdir, edir=str(edir), idir=idir )
            ipop = {}
            
            with open( myfile, 'r' ) as infile:

                key_heading = infile.readline()
                keys   = key_heading.split()
                values = [ [] for ivalue in range(len(keys)) ]
                ipop = dict( zip(keys,values) )
                
                all_lines = infile.readlines()
                for iline in all_lines :
                    ilist = iline.split()
                    for i in range(len(keys)) : ipop[keys[i]].append( float(ilist[i]) )
            
            self.pop[idir] = copy.deepcopy( ipop )

             

    def fit_func( self, angles, *params ) :

        theta, phi = angles

        def ftheta( angle ) :
            return [ np.sin(angle)**2, 
                     np.sin(angle)**4,
                     np.cos(angle) * np.sin(angle)**2, 
                     np.cos(angle) * np.sin(angle)**4, 
                     np.sin(angle)**6 ]            

        def fphi( angle ) :
            fphi_list  = [ np.cos(i*angle) for i in range(1,6) ]
            fphi_list += [ np.sin(i*angle) for i in range(1,6) ]
            fphi_list += [ np.cos(6*angle) ] 
            return fphi_list

        ftheta_list = ftheta( theta )
        fphi_list   = fphi( phi )

        fvalue = 0.0
        for i in range(7) : fvalue = fvalue + params[i] * np.cos(theta)**i
            
        ip = 6
        for itheta in ftheta_list :
            for iphi in fphi_list :
                ip   += 1
                fvalue = fvalue + params[ip] * itheta * iphi        

        return fvalue



    def get_angular( self, itime=0 ) :

        import scipy.optimize as optimization
    
        #: theta  and phi list
        theta = [ self.field[idir]['theta0']*deg2rad for idir in range(1,self.ndir+1) ]
        phi   = [ self.field[idir]['phi0']*deg2rad   for idir in range(1,self.ndir+1) ]
        
        #: linear parameters to optimize, initialize to one 
        params = np.ones( self.ndir )

        #: data to fit
        if itime == 0 : itime = self.ntimes - 1
        ref_data = [ 1.0 - self.results[idir]['norm2'][itime] for idir in range(1,self.ndir+1) ]
        
        popt, pcov = optimization.curve_fit( self.fit_func, (theta,phi), ref_data, params, method='lm' )
        
        #: print out fitted values 
        i = 0
        for (itheta,iphi) in zip(theta,phi) :
            myval = self.fit_func(  (itheta,iphi), *popt )
            print( "ref_val = {:15.10f}    fitted_val = {:15.10f}".format( ref_data[i], myval ) )
            i+=1
        return popt
        
        
        
    def plot_angular(self, popt, getxyz = False ) :

        #from mpl_toolkits.mplot3d import Axes3D
        #from matplotlib.colors import LightSource

        #: set grid
        theta, phi = np.linspace(0, np.pi, 40), np.linspace(0, 2.0*np.pi, 40)

        #: get data on mesh grid
        fitted_data = np.zeros( (len(theta),len(phi)) )
        for itheta in range(theta.size) :
            for iphi in range(phi.size) :
                angles = ( theta[itheta], phi[iphi] )
                fitted_data[iphi,itheta] = np.absolute( self.fit_func( angles, *popt ) )
                
        #: get points to plot
        theta_mesh, phi_mesh = np.meshgrid( theta, phi )
        x = np.sqrt(fitted_data) * np.cos(phi_mesh) * np.sin(theta_mesh)
        y = np.sqrt(fitted_data) * np.sin(phi_mesh) * np.sin(theta_mesh)
        z = np.sqrt(fitted_data) * np.cos(theta_mesh)

        if getxyz : return x, y, z
        
        #: lighting source
        light = LightSource(270, 90)
        rgb = np.ones((z.shape[0], z.shape[1], 3))
        pink = [ 219/255, 112/255, 147/255 ]
        illuminated_surface = light.shade_rgb(rgb*pink, z)

        #: plot plot plot plot
        fig = plt.figure()
        ax  = plt.axes( projection='3d' )
        ax.plot_surface( x,y,z,alpha=0.4,linewidth=0,antialiased=False,facecolors=illuminated_surface)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        #ax.set_axis_off()
        plt.show()


    def plot_contour(self, popt ) :
    
        x = {'30':[], '60':[], '90':[] }
        y = {'30':[], '60':[], '90':[] }
        z = {'30':[], '60':[], '90':[] }
        phi = np.linspace(0, 2.0*np.pi, 40)
        for itheta in ( 30, 60, 90 ) :
            theta = itheta*deg2rad
            fitted_data = np.zeros( (40 ) )
            for iphi in range(phi.size) :
                angles = (theta,phi[iphi])
                fitted_data[iphi] = np.absolute( self.fit_func( angles, *popt ) )
                x[str(itheta)] = np.sqrt(fitted_data) * np.cos( phi ) #* np.sin(theta)
                y[str(itheta)] = np.sqrt(fitted_data) * np.sin( phi ) #* np.sin(theta)
                z[str(itheta)] = np.sqrt(fitted_data) * np.cos(theta)
                
        fig, ax1 = plt.subplots()
        ax1.plot( x['30'], y['30'], label='theta=30' )
        ax1.plot( x['60'], y['60'], label='theta=60' )
        ax1.plot( x['90'], y['90'], label='theta=90' )
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.legend()
        plt.show()


#---------------------------------------------------------------------------------------------------#
#ch3i = Results()
#ch3i.ndir = 62
#ch3i.ntimes = 320

#ch3i.read_results( rdir='/wsu/home/gq/gq39/gq3921/SYSTEMS/ch3i/static/RESULTS/',edir='1' )
#popt = ch3i.get_angular(itime=200)
#ch3i.plot_contour(popt)
