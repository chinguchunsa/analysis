#!/wsu/el7/pre-compiled/python/3.7/bin/python

import numpy  as np
import orbkit as ok
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import animation


#plt.switch_backend('agg')
mpl.rcParams['axes.linewidth'] = 2.0     #: thicker border
mpl.rcParams['xtick.major.width'] = 2.0  #: thicker xtic marks
mpl.rcParams['ytick.major.width'] = 2.0  #: thicker ytic marks


myproc = 16


class TDCI_Wavefunction( ok.QCinfo ) :

    def __init__( self, f_chkpt='', f_moden='' ) :

        """ Initialize variables.  Inherit QCinfo """

        ok.QCinfo.__init__(self)
        self.ntimes = 0         #: total number of timesteps
        self.orbs   = {}        #: orbital info used for TDCI (may not include core)
        self.f_chkpt = f_chkpt  #: checkpoint file REQUIRED
        self.f_moden = f_moden  #: mo_den file REQUIRED.  Generated with analysis.exe
        self.time   = []        #: list of timepoints
        self.norm2  = []        #: list of norms
        self.moden_info = []    #: list of dictionary { 'mo_pair':(i,j), 'coeff2': float }
        self.mo_ongrid = []     #: list of lists.  each list contains list of dictionaries
        self.rho0      = []     #: rho0.  needed to get difference density.  must specify beforehand to get difference density!
        self.contour   = {}     #: contour variables
                


    def set_standard_grid(self) :

        ok.grid.min_ = [-8.0, -8.0, -8.0]   #: Specifies minimum grid values (regular grid).
        ok.grid.max_ = [ 8.0,  8.0,  8.0]   #: Specifies maximum grid values (regular grid).
        ok.grid.N_   = [ 101,  101,  101]   #: Specifies the number of grid points (regular grid).
        ok.grid.grid_init(is_vector=False, force=False)

        
    def set_logfile(self, f_log='none' ) :
        
        ok.options.outputname = f_log

    
    def initialize( self, grid=True ) :

        """ 'Officially' initialize class and assign class variables from checkpoint file 
                        There must be a more python-ic way to do this.  
        """
        
        qc_tmp = ok.main_read( self.f_chkpt, itype='gaussian.fchk', all_mo=True )

        self.geo_info     = qc_tmp.geo_info
        self.geo_spec     = qc_tmp.geo_spec
        self.ao_spec      = qc_tmp.ao_spec
        self.ao_spherical = qc_tmp.ao_spherical
        self.mo_spec      = qc_tmp.mo_spec
        self.etot         = qc_tmp.etot
        self.com          = qc_tmp.com
        self.coc          = qc_tmp.coc
        self.bc           = qc_tmp.bc
        self.pop_ana      = qc_tmp.pop_ana
        self.states       = qc_tmp.states
        self.alpha_homo   = qc_tmp.alpha_homo
        self.beta_homo    = qc_tmp.beta_homo
        self.dipole_moments = qc_tmp.dipole_moments        
        del qc_tmp

        self.set_logfile( f_log = 'OK_log' )
        if grid : self.set_standard_grid()

        

    def write_vmd_tcl( self, start, end, outfile, diff, isoval=0.0005 ) :
        
        """ Write VMD script to animate movie """

        import datetime as dt

        repids   = [ '$repid1' ]
        updreps  = [ 'updrep1' ]
        isovals  = [ isoval ]
        colorids = [ 0 ]               #blue
        if diff : 
            repids.append( '$repid2')
            updreps.append( 'updrep2' )
            isovals.append( -isoval )  #loss
            colorids.append( 1 )       #red      

        f_vmd = open( 'VMD_' + outfile[:-2]+'.tcl', 'w' )
        
        mystr = '# generated @ ' + str(dt.datetime.now()) + '\n'
        mystr = mystr + '# source me in VMD terminal\n\n'
        #: files to be loaded 
        mystr = mystr + 'set updmol [mol new {{{0}}} type cube waitfor all]\n'.format( outfile+str(start)+'.cube' )
        mystr = mystr + 'for {{ set i {0} }} {{ $i <= {1} }} {{ incr i }} {{\n'.format( start+1, end )
        mystr = mystr + '\t set myfile {0} \n'.format(outfile+'$i'+'.cube')
        mystr = mystr + '\t mol addfile $myfile type cube waitfor all \n' 
        mystr = mystr + '}\n\n' 
        #: define how to draw molcule 
        mystr = mystr + 'mol delrep 0 top\n' 
        mystr = mystr + 'mol representation CPK 1.700000 0.300000 20.000000 16.000000\n' 
        mystr = mystr + 'mol color Name\n' 
        mystr = mystr + 'mol addrep top\n\n' 
        #: isosurface 
        for i in range(len(updreps)) :
            mystr = mystr + 'mol representation Isosurface {0} 0.0 0.0 0.0\n'.format( isovals[i] ) 
            mystr = mystr + 'mol Material Transparent\n' 
            mystr = mystr + 'mol color ColorID {0}\n'.format( colorids[i] ) 
            mystr = mystr + 'mol addrep top\n\n' 
        #: update 
        for i in range(len(updreps)) : mystr = mystr + 'set {0} [mol repname top {1}]\n'.format( updreps[i], i+1 ) 
        #: update movie
        mystr = mystr + '\nproc update_iso {args} {\n' 
        mystr = mystr + '   global updmol\n' 
        for i in range(len(updreps)) : mystr = mystr + '   global {0}\n'.format( updreps[i] ) 
        mystr = mystr + '\n' 
        for i in range(len(updreps)) :
            mystr = mystr + '   set {0} [mol repindex $updmol {1}]\n'.format( repids[i][1:],'$'+updreps[i] ) 
            mystr = mystr + '   if {{ {0} < 0}} {{ return }}\n'.format( repids[i] ) 
        mystr = mystr + '\n   set frame [molinfo $updmol get frame]\n\n' 
        for i in range(len(updreps)) :        
            mystr = mystr + '   lassign [molinfo $updmol get "{{rep {0}}}"] rep\n'.format( repids[i] ) 
            mystr = mystr + '   mol representation [lreplace $rep 2 2 $frame]\n' 
            mystr = mystr + '   mol color ColorID {0}\n'.format( colorids[i] )
            mystr = mystr + '   mol modrep {0} $updmol\n\n'.format( repids[i] )
        mystr = mystr + '}\n\n' 
        mystr = mystr + 'trace variable vmd_frame($updmol) w update_iso\n' 
        mystr = mystr + 'animate goto 0\n' 

        f_vmd.write( mystr )
        f_vmd.close()



    def read_f_moden( self ) :
        
        """ Reads MO-e*-d*.out file produced from ANALYSIS.f90 
            NOTE  orbital indexing starts from 0
        """

        import copy as cp

        #: initialize 
        self.moden_info = []

        with open( self.f_moden, 'r' ) as infile : 
        
            #: read first four lines 
            self.ntimes = int( infile.readline().split()[-1] )
            self.orbs['alpha_homo'] = int( infile.readline().split()[-1] )
            self.orbs['beta_homo']  = int( infile.readline().split()[-1] )
            self.orbs['alpha_tot']  = int( infile.readline().split()[-1] ) 
            
            #: compute number of core orbitals not included in TDCI
            self.orbs['alpha_offset'] = self.alpha_homo - self.orbs['alpha_homo']
            self.orbs['beta_offset']  = self.beta_homo  - self.orbs['beta_homo']

            #: shorter named variables for pretty code
            dalpha = self.orbs['alpha_offset']
            dbeta  = self.orbs['alpha_offset'] + self.orbs['beta_offset']             

            
            #: initialize. set key 
            iline, self.time, self.moden_info, self.norm2 = 0, [], [], []
            tmpdict, key  = {}, ( 'mo_pair', 'coeff2' )
            

            #: read in rest of the file
            all_lines = infile.readlines() 

            #: parse lines
            while iline < len(all_lines)  :
                
                #: initialize temp array
                itime_moden_info = []
                
                self.time.append( float(iline) )

                #: include core
                if self.orbs['alpha_offset'] != 0 :
                    for i in range( self.orbs['alpha_offset'] ) :
                        mo_pair, coeff2 = (i,i) , 1.0
                        tmpdict = dict(zip( key, (mo_pair,coeff2) ))
                        itime_moden_info.append( cp.deepcopy(tmpdict) )
                    for i in range( self.orbs['beta_offset'] ) :
                        ii = self.orbs['alpha_offset'] + self.orbs['alpha_tot'] + i 
                        mo_pair, coeff2 = (ii,ii), 1.0
                        tmpdict = dict(zip( key, (mo_pair,coeff2) ))
                        itime_moden_info.append( cp.deepcopy(tmpdict) )

                iline += 1

                #: parse list 
                curr_line = all_lines[iline]
                for mygroup in curr_line.split(',') : 
                    try :
                        my_mo1, my_mo2, my_coeff = mygroup.split()
                        #: -1 cause python starts from 0, fortran from 1
                        my_mo1 , my_mo2 = int(my_mo1)-1, int(my_mo2)-1
                        #: add core offset
                        my_mo1  = my_mo1 + ( dalpha if my_mo1 < self.orbs['alpha_tot'] else dbeta )
                        my_mo2  = my_mo2 + ( dalpha if my_mo2 < self.orbs['alpha_tot'] else dbeta )
                        mo_pair, coeff2 = ( my_mo1 , my_mo2 ), float(my_coeff)
                        tmpdict = dict(zip( key, (mo_pair,coeff2) ))
                        itime_moden_info.append( cp.deepcopy(tmpdict) )
                    except : pass                    
                        
                iline += 1
                
                #: list_of_dictionaries_@itime = moden_info[itime]
                self.moden_info.append( itime_moden_info )

                

    def get_density_ground ( self ) :

        """ Assign rho ground """

        self.rho0 = ok.rho_compute( self, itype='gaussian.fchk', calc_mo=False, numproc=myproc )



    def get_itime_density( self, itime=0, diff=False ) :

        """ Generate density at time itime """

        itime_moden_info = self.moden_info[itime]
        nterms  = len( itime_moden_info )

        (i0,j0) = itime_moden_info[0]['mo_pair']
        i0_mo   = self.mo_ongrid[i0,:,:,:]

        gridx, gridy, gridz = i0_mo.shape
        density = np.zeros( (gridx,gridy,gridz) )        

        #: go through mo list
        for idensity in range( nterms ) :
            ( i,j )  = itime_moden_info[idensity]['mo_pair']
            my_coeff = itime_moden_info[idensity]['coeff2']
            if i != i0 :
                i0 = i
                i0_mo = self.mo_ongrid[i0,:,:,:]
            density += my_coeff * i0_mo * self.mo_ongrid[j,:,:,:]

        return ( density-self.rho0 ) if diff else density



    def generate_density_movie( self, start=-100, end=-100, outfile='default', diff=False, write_vmd=True ) :
        
        """ Generate cube files 
            or send back density to generate_contour
        """

        #: set outfile prefix
        if outfile == 'default' : 
            outfile = self.f_moden.replace('MO','DEN') if 'MO' in self.f_moden else self.f_moden
            if '.out' in outfile : outfile = outfile[:-4] 

        #: read MO indices and product of MO coefficients 
        self.read_f_moden()

        #: get all MOs and rho0
        self.mo_ongrid = ok.rho_compute( self, itype='gaussian.fchk', calc_mo=True, numproc=myproc )
        
        # get start time and end time
        if start == -100 : start = 0
        if end   == -100 : end   = len( self.time ) - 1        
            
        # write vmd TCL file
        isoval = ( 0.0005 if diff else 0.001 )
        ofile = outfile + str(start)
        if write_vmd : self.write_vmd_tcl( start, end, ofile, diff, isoval )

        for itime in range( start, end+1 ) :                
            ofile   = outfile + '-t' + str(itime)
            density = self.get_itime_density( itime=itime, diff=diff )
            ok.main_output( density, self.geo_info, self.geo_spec, outputname=ofile, otype='cb' )



    def read_cube( self, f_cube='' ) :

        """ Read in cube files 
        Returns density on grid
        """

        with open( f_cube,'r' ) as infile :
            
            #: read title
            for i in range(2) : infile.readline()
            #: read in min vals
            minline = [ ival for ival in infile.readline().split() ]
            xmin, ymin, zmin = float(minline[1]), float(minline[2]), float(minline[3])
            #: read x grid points
            xline = [ ival for ival in infile.readline().split() ]
            nx_pts, dx = int(xline[0]), float(xline[1])
            #: read y grid points
            yline = [ ival for ival in infile.readline().split() ]
            ny_pts, dy = int(yline[0]), float(yline[2])
            #: read z grid points
            zline = [ ival for ival in infile.readline().split() ]
            nz_pts, dz = int(zline[0]), float(zline[3])

            
            #ok.grid.min_ = [ xmin, ymin, zmin ]
            #ok.grid.delta_ = [ dx, dy, dz ]
            #ok.grid.max_ = [ xmin + dx * nx_pts, 
            #                 ymin + dy * ny_pts, 
            #                 zmin + dz * nz_pts ]
            #ok.grid.N_ = [ nx_pts, ny_pts, nz_pts ]            
            #ok.grid.grid_init(is_vector=False, force=False)
            print(ok.grid.get_grid())
            
            
            for i in range( len(self.geo_info) ) : infile.readline()

            density = []
            for iline in infile.readlines() :
                for ielement in iline.strip().split() : density.append( float(ielement) )
            density = np.array( density )
            
        return density.reshape( (len(ok.grid.x),len(ok.grid.y),len(ok.grid.z)) )


    
    def get_diff_density_from_cube( self, start=-100, end=-100, cfiles='none', ofiles='none' ) :

        """ Read in cube files 
        Returns difference density cube files 
        Assumes standard grid
        Must set self.rho0 externally
        Must provide cfiles = [ 'file1.cube', 'file2.cube', .. ]
        """

        try :
            if self.rho0[0] == 'NA' : 
                print( ' ERROR:  must specify rho0 first' )
                print( ' for example with               ' )
                print( '     obj.rho0 = obj.read_cube( "cubefile.cube" )' )
                return
        except : pass

        #: get cube files
        if cfiles == 'none' : 
            print( ' ERROR:  must specify list of cube files' )
            print( ' for example with one file, cfiles=["DEN-e1-d1-t1.cube"]' )
            return

        #: outfiles prefixes or list
        if ofiles == 'none' : 
            print( ' ERROR:  must specify ofiles = [ output_cubefiles_name ] ' )
            print( ' or ofiles = output_cubefile_prefix ' )
            return
                   
        # get start time and end time
        if start == -100 : start = 0
        if end   == -100 : end   = len( self.time ) - 1

        for i in range( len(cfiles ) ) :
            rho = self.read_cube( f_cube=cfiles[i] )
            if type(ofiles) == list : 
                output = ofiles[i][:-5] if '.cube' in ofiles[i] else ofiles[i]
            else :
                output = ofiles + '-t' + str(i)    
            ok.main_output( rho - self.rho0, self.geo_info, self.geo_spec, outputname=output, otype='cb' )
            


    def setup_contour( self, xplane=False, yplane=False, zplane=False ) :
        
        """ Set up contour variables """

        self.contour['min1'] = -8
        self.contour['max1'] = 8
        self.contour['min2'] = -8 
        self.contour['max2'] = 8

        self.contour['xgrid'] = ok.grid.x
        self.contour['xmin']  = ok.grid.min_[0]
        self.contour['xmax']  = ok.grid.max_[0]
        
        self.contour['ygrid'] = ok.grid.y
        self.contour['ymin']  = ok.grid.min_[1]
        self.contour['ymax']  = ok.grid.max_[1]

        self.contour['zgrid'] = ok.grid.z
        self.contour['zmin']  = ok.grid.min_[2]
        self.contour['zmax']  = ok.grid.max_[2]

        self.contour['dx'] = (self.contour['xmax'] - self.contour['xmin']) / float( len(self.contour['xgrid']) )
        self.contour['dy'] = (self.contour['ymax'] - self.contour['ymin']) / float( len(self.contour['ygrid']) )
        self.contour['dz'] = (self.contour['zmax'] - self.contour['zmin']) / float( len(self.contour['zgrid']) )
            
        #:indices
        self.contour['i_x0'] = len(self.contour['xgrid'])//2  #: index for x=0
        self.contour['i_y0'] = len(self.contour['ygrid'])//2  #: index for y=0
        self.contour['i_z0'] = len(self.contour['zgrid'])//2  #: index for z=0
        
        #:indices
        if xplane : min1, min2, max1, max2, dxyz1, dxyz2 = 'zmin', 'ymin', 'zmax', 'ymax', 'dz', 'dy'
        if yplane : min1, min2, max1, max2, dxyz1, dxyz2 = 'zmin', 'xmin', 'zmax', 'xmax', 'dz', 'dx'
        if zplane : min1, min2, max1, max2, dxyz1, dxyz2 = 'ymin', 'xmin', 'ymax', 'xmax', 'dy', 'dx'
        self.contour['i_min1'] = int( ( self.contour['min1'] - self.contour[min1] ) / self.contour[dxyz1] )
        self.contour['i_max1'] = int( ( self.contour['max1'] - self.contour[min1] ) / self.contour[dxyz1] )
        self.contour['i_min2'] = int( ( self.contour['min2'] - self.contour[min2] ) / self.contour[dxyz2] )
        self.contour['i_max2'] = int( ( self.contour['max2'] - self.contour[min2] ) / self.contour[dxyz2] )
        

    def generate_contour( self, start=-100, end=-100, diff=False, xplane=False, yplane=False, zplane=False, outfile='default.gif' ) :
        
        """ Generate contour plots along specified plane 
        """
        
        #: get "CI/MO" info
        self.read_f_moden()

        #: get MO on grid
        self.mo_ongrid = ok.rho_compute( self, itype='gaussian.fchk', calc_mo=True, numproc=myproc )
        self.get_density0()

        #: set up contour
        self.setup_contour( xplane, yplane, zplane )
        
        #: pretty variables
        imin1, imin2 = self.contour['i_min1'], self.contour['i_min2']
        imax1, imax2 = self.contour['i_max1'], self.contour['i_max2']
        ix0, iy0, iz0 = self.contour['i_x0'], self.contour['i_y0'], self.contour['i_z0']
        
        # get start time and end time
        if start == -100 : start = 0
        if end   == -100 : end   = len( self.time ) - 1
        frames = [ i for i in range( start, end + 1 ) ]

        #: initialize 
        fig, ax  = plt.subplots()
        mylevels = ( np.arange(-0.01,0.01,0.00001) if diff else np.arange(0,0.03,0.00001) )       

        if xplane :
            grid1 = self.contour['zgrid'][ imin1:imax1 ]
            grid2 = self.contour['ygrid'][ imin2:imax2 ]
            density = self.get_itime_density(0,diff=diff)[ ix0, imin2:imax2, imin1:imax1 ] 
            ax.set_xlabel('z') 
            ax.set_ylabel('y')
        if yplane :
            grid1 = self.contour['zgrid'][ imin1:imax1 ]
            grid2 = self.contour['xgrid'][ imin2:imax2 ]
            density = self.get_itime_density(0,diff=diff)[ imin2:imax2, iy0, imin1:imax1 ]
            ax.set_xlabel('z') 
            ax.set_ylabel('x')
        if zplane :
            grid1 = self.contour['ygrid'][ imin1:imax1 ]
            grid2 = self.contour['xgrid'][ imin2:imax2 ]
            density = self.get_itime_density(0,diff=diff)[ imin2:imax2, imin1:imax1, iz0 ] 
            ax.set_xlabel('y') 
            ax.set_ylabel('x')
        init = ax.contourf( grid1, grid2, density, mylevels, cmap='RdGy' )
        fig.colorbar(init)

        #: animate
        def animate(itime) :
            ax.clear()
            if xplane : 
                density = self.get_itime_density(itime,diff=diff)[ ix0, imin2:imax2, imin1:imax1 ]
                ax.set_xlabel('z')
                ax.set_ylabel('y')
            if yplane : 
                density = self.get_itime_density(itime,diff=diff)[ imin2:imax2, iy0, imin1:imax1 ]
                ax.set_xlabel('z')
                ax.set_ylabel('x')
            if zplane : 
                density = self.get_itime_density(itime,diff=diff)[ imin2:imax2, imin1:imax1, iz0 ]
                ax.set_xlabel('y')
                ax.set_ylabel('z')
            ax.contourf( grid1, grid2, density, mylevels, cmap='RdGy' )
            ax.set_title(itime)

        ani = animation.FuncAnimation( fig, animate, frames=frames, interval=50, blit=False, repeat=False )
        ani.save( outfile )  
        

    def get_mo_ongrid( self, mo_id='none', nproc=myproc ) :
        
        if type(mo_id) != int :
            print( ' ERROR:  must provide MO # ' )
            return

        self_dict = self.todict()
        self_dict['mo_spec'] = []
        self_dict['mo_spec'].append( self.mo_spec[mo_id] )

        return ok.rho_compute( self_dict, calc_mo=True, numproc=nproc )
    

    
    def write_cube( self, rho=0, cfile='default' ) :
        
        if '.cube' in cfile : cfile = cfile[:-5]
        ok.main_output( rho, self.geo_info, self.geo_spec, outputname=cfile, otype='cb' )



#:-------------------------------------------------------------------------------------------------------------------------:#
#: molecule = TDCI_Wavefunction( f_moden='file_name', f_chkpt='chkpoint.fchk' )     REQUIRED
#: molecule.set_logfile( f_log = 'file_name' )
#: molecule.initialize( grid=True )                                                 REQUIRED
#: molecule.write_vmd_tcl( start=start_time, end=end_time, outfile='script_name', diff=True_False, isoval=0.001 )
#: molecule.read_f_moden( normalize=True_False )
#: molecule.get_density_ground()                                                    RETURNS rho[ngrid_x,ngrid_y,ngrid_z]
#: molecule.get_itime_density( itime=itime, diff=True_False )                       RETURNS rho or rho - molecule.rho0
#: molecule.generate_density_movie(start=start_time, end=end_time, outfile='prefix_for_cubefiles', normalize=True_False, write_vmd=True_False )
#: molecule.read_cube( f_cube='cubefile.cube' )                                     RETURNS rho[ngrid_x,ngrid_y,ngrid_z]
#: molecule.get_diff_density_from_cube( start=start_time, end=end_time, cfiles=[LIST OF CUBEFILES], ofiles=[LIST OF OUTPUT]_or_prefix )
#: molecule.get_mo_ongrid( mo_id=#, nproc=numproc )                                RETURNS rho
#: molecule.write_cube( rho = something_on_grid, cfile='cube_file_name' )
#:-------------------------------------------------------------------------------------------------------------------------:#


n2 = TDCI_Wavefunction( f_moden='../DENSITY_MOVIE-e12-d1.out', f_chkpt='/wsu/home/gq/gq39/gq3921/SYSTEMS/n2/gaussian/cis/RA_6_5/n2.fchk' )
n2.initialize( grid=True )
n2.get_density_ground()
n2.generate_density_movie( outfile='diff-e12-d1', diff=True  )
