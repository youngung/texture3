"""
This module started as a hobby project while learning Python back in 2011
when I was a graduate student. At the moment I am using it for my custom pole
figures from discrete orientations used for VPSC-like crystal plasticity
codes inclduing VPSC, dEVPSC

Features:
  crystal symmetries : cubic and hexagonal
  This module can be easily extened for other crystal symmetries.
  At least, it additionally works for orthorhombic crystal structure
  for a paper I worked on uranium. The crystal symmetry should be
  extended and generalized for all crystal symmetries existing, which
  has not been pursued yet.


Pole figure and inverse pole figure plotting by stereographic
projection. Both contour and dot types are available for pole
figures, whereas only dot type is available for inverse pole figure, yet.

It should be machine independent. Make the way treating the data
uniform and intiutive and simple.

--------
examples
--------
 >>> import upf
 >>> mypf = upf.polefigure(ngrain=8000,
                           # or filename ='texture/08000.cmb'
                           csym = 'cubic', cdim=[1.,1.,1.],
                           cang =[90.,90.,90.])
 >>> cnts = mypf.pf(pole=[ [1,0,0],[1,1,0],[1,1,1]],
                    mode='contourf', ifig=2,
                    dm=7.5, dn=7.5, levels=None, cmode='gray_r')

        * Pole figure plotting upon the polycrystal aggregate.
          --> refer to 'class polefigure'

        * Experimental pole figure plotting software where
         the range of khi (tiliting angle) is often limited.
         However, the plotting scheme will be very similiar to the
         pole figure for polycrystalline aggregates, in that
         experimental pole figure often comes as a grid set of data
         in the polar coordinate system (theta as angle, and khi as
         radius)
          --> First of all, I should come up with a uniform format
          having a grid format consisting of theta and khi axes.
          the resolution of grid should be decided after getting the
          data.

          range of data:
          theta : 0. ~ 2.*pi
          khi : 0. ~ pi/2.

          detecting the np.shape and decide the incremental resolution
          of each axis.

          Let's say it will be something like below:

          >>> m, n = np.shape(data)
               # m is the number of grid along theta axis
               # while n is the number of grid along khi

          --> Often, as far as I know, experimental pole figure should be
          normalized since the incremental scanning area is varying with
          increasing khi(tilting angle). The issues related to background
          and defocus effect should be dealt with independently, as well.

          The format controller is also necessary.
          Current cases of interests:
              1. Bruker D8 Discovery (at GIFT POSTECH) (UXD)
              2. Unknown system at Dr. Steglich's institute
              3. Experimental PF (popLA format)

 >>> import upf
 >>> mypf = upf.polefigure(epf='???.epf')
..
..
 >>> cnts = mypf.pf()
"""
# print __doc__

"""
Updates logs are no longer tracked as the project is version-controlled using Git.
  # 1
  2011-16-Mar
  Contour plot is done on the given axes if necessary.

  # 2
  2011-17-April
  The symmetry module is cythonized.
  The __equiv__ method has been taken out from the loop over grains,
  which saves some computational waste.

  # 3
  2011-18-April
  The hardwired smoothing rim is expanded to the 2nd rim.
  Thus now it is n=0, n=1 axes along which the intensities are averaged.

  # 4
  2011-14-Sept
  The pixel view superimposition to contour pole figure is initiated.
  Cython part having been used for 'sym.py' has been removed.
"""
#----------------------------------------------------------------------c
## import library blocks
import warnings
# warnings.filterwarnings("ignore")
# warnings.filterwarnings("error")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib #matplotlib as raw
import os, glob, math
from .randomEuler import randomEuler as re
from .euler import euler # in euler module def euler:
                        # A-matrix and Euler angles
import time,random
import fortranformat as ff
import sys

try:
    import MP
except:
    print('-'*50)
    print('MP was not installed')
    print('You may clone it and install via:')
    print('git@github.com:youngung/mpl-lib.git')
    print('-'*50)
else:
    from MP import progress_bar
    t2s = progress_bar.convert_sec_to_string
    uet = progress_bar.update_elapsed_time

## removing joblib
is_joblib=False
pi   = math.pi;
cos  = math.cos; sin  = math.sin

# import pp ## parallel

def cub(filename=None,gr=None,ifig=3,**kwargs):
    """
    arguments
    =========
    filename = None
    gr       = None
    ifig     = 3
    pole     = [[1,0,0],[1,1,0],[1,1,1]],
    cmode    = None):
    """
    if filename!=None: mypf = polefigure(filename=filename, csym='cubic')
    elif gr!=None: mypf = polefigure(grains=gr, csym='cubic')
    fig=mypf.pf_new(**kwargs)
    return fig

def cubgr(gr=None,ifig=3,poles=[[1,0,0],[1,1,0],[1,1,1]]):
    """
    Arguments
    =========
    gr
    ifig
    poles
    """
    mypf = polefigure(grains=gr, csym='cubic')
    fig = mypf.pf_new(poles=poles,cmap='jet')
    return fig

def pfnorm(data):
    """
    experimental incomplete pole figure preliminary normalization

    data format should be correct
    grid (phi, khi)

    resolution should be 5.0 degress in both phi and khi
    phi range: 0~ 355
    khi range: 0~ 5*(nn-1)

    Argument
    ========
    data
    """
    # All angles are in radian
    if len(data)!=72:
        print('number of phi grid:  %i'%len(data))
        raise IOError('Unexpected resolution along phi axis')

    dphi = 360. / len(data)
    dphi = dphi * np.pi / 180. # dphi
    dkhi = 5. * np.pi / 180. # dkhi
    print('dkhi, dphi', dphi*180./np.pi, dkhi*180./np.pi)

    nkhi = len(data[0])
    phi_i = 0.
    phi_f = np.pi * 2.
    khi_i = 0.
    khi_f = dkhi * (nkhi - 1)
    print('khi range', khi_i, khi_f*180./np.pi)
    print('phi range', phi_i, phi_f*180./np.pi)
    # spanned area, i.e., the area of incomplete hemisphere:
    # area = (np.cos(khi_f) - np.cos(khi_i)) * (phi_f - phi_i)

    Nf = 0.

    a = 0.
    b = 0.
    ## below needs vectorization
    for i in range(len(data)):
        for j in range(len(data[0])):
            a = a + np.sin(dkhi*j)

    for i in range(len(data)):
        for j in range(len(data[i])):
            b = b + data[i,j] * np.sin(dkhi*j)

    # for k in xrange(len(data)):
    #     for l in xrange(len(data[k])):
    #         b = b + data[k, l] * np.sin(dkhi*l)

    for i in range(len(data)): #phi
        for j in range(len(data[i])): #khi
            data[i,j] = data[i,j] * a / b

    return data

def epfformat(mode=None, filename=None):
    """
    Experimental pole figure format controller
    mode:
      "steglich"
      "bruker"  *.uxd file
      "epf" (2011-Oct-6) epf popLA experimental pole figure format
      "xpc" (2017-Feb) xcp format compliant with MAUD

    Returns the pole figure data as
    the standard format (m x n numpy array)
    each of axes stands for rotating (phi) and tilting (khi)
    angle in the laboratory space.
    These angles will be angle and radius in the
    space onto which a pole figure is projected.

    conventions:
     tilting angle: khi
     rotating angle: phi
     dk: incremental khi angle
     dp: incremental phi angle
     nk: number of points along khi axis
     np: number of points along phi axis

     angles are converted into radian whenever possible

    Arguments
    =========
    mode     = None
    filename = None

    Returns
    =======
    data
    max_khi
    hkl
    """

    ## steglich's format
    nskip_calc = 13
    nskip_raw = 28
    ##

    if mode=='steglich':
        # Calculated or raw pole figure
        i = 0
        print('filename=', filename)
        if os.path.isfile(filename)!=True:
            raise IOError('file is not available')

        while True:
            try:
                data = np.loadtxt(
                    filename, skiprows=i)
            except: i = i + 1
            else:
                print('number of skipped rows: %i'%i)
                break
            if i>1000: raise IOError('something is wrong')

        ## raw pole figure format
        if i==nskip_raw:
            # axes: (phi, khi)
            data = data.T # (khi, phi)

            # upon each axis
            f = open(filename)
            temp = f.readlines()
            khi = list(map(float, temp[nskip_raw - 1].split()[1:]))
            khi = np.array(khi)
            khi = khi * np.pi/ 180.
            phi = data[0] #first column is for phi
            phi = phi * np.pi/ 180.
            data = data[1:] #remove phi
            data = data.T # (phi, khi) axis

            ## dp and dk
            dp = phi[1] - phi[0]
            dk = khi[1] - khi[0]

            ## shift phi back-modification
            phi, data = shiftphi(data, dp, phi)
            ## shift phi back-modification ends

            ## check if there is phi=0, 360 at the same time
            ## if that's the case remove data along phi=360
            phi, data = duplicate0_360(phi, data)
            ##

            ## check if this is incomplete pole figure
            tiny = 0.000001
            isincomplete=False
            if np.pi/2. - khi[-1] > tiny:
                print('Incomplete pole figure')
                print('khi range: %3.2f ~%3.2f'%(
                    khi[0]*180./np.pi, khi[-1]*180./np.pi))
                isincomplete=True

                ## normalization
                if input('y(norm), or n(no)>>>')=='y':
                    data = pfnorm(data)

                dum = np.zeros((np.pi*2./dp, np.pi/2./dk+1))
                for i in range(len(data)):
                    for j in range(len(data[i])):
                        dum[i, j] = data[i, j]
                data = dum.copy()
                del dum

        ## calculated pole figure format
        elif i==nskip_calc: #He had two formats: raw and calculated.
            print('Calculated pole figure format')
            # axes: (phi, khi)
            data = data.T #(khi, phi)
            f = open(filename)
            temp = f.readlines()
            khi = list(map(float, temp[nskip_calc - 1].split()[1:]))
            khi = np.array(khi)
            khi = khi * np.pi / 180.
            phi = data[0]
            phi = phi * np.pi / 180.
            data = data[1:]
            data = data.T #(phi, khi)

            ## dp and dk
            dp = phi[1] - phi[0]
            dk = khi[1] - khi[0]

            ## shift phi back-modification
            phi, data = shiftphi(data, dp, phi)
            ##

            ## check if there is phi=0, 360 at the same time
            phi, data = duplicate0_360(phi, data)
            ##

    elif mode=='bruker':
        ## make full use of existing uxd.py script
        ## --> normalization is missing...
        ## This must be completed!
        from . import uxd
        print('You are now in the bruker mode under epfformat')
        print('given file name is %s'%filename)
        myuxd = uxd.pf(filename=filename, mode='pf')
        if len(myuxd.polefigures)>1:
            print('multiple pole figures are found')
            input()
        for i in range(len(myuxd.polefigures)):
            pf = myuxd.polefigures[i]
            pf = pfnorm(pf) ##normalize

    elif mode=='epf':
        """
        ready made popLA epf format parser
        consider the possibility of multiple number of polefigure
        ## phi must be 0~355, khi must be 0~90 with 5 as an ang
        resolution
        """
        print('You are now reading experimental pole figure(s) :%s'%filename)
        blocks = open(filename, 'rU').read().split('(')[1:]
        if len(blocks)==0:
            msg1 = 'epf parser in upf assumes that hkl is embraced by paratheses'
            msg2 = ' %s has no paranthesis that embraces hkl was found'%filename
            msg  = '%s \n %s'%(msg1,msg2)
            print('-'*52)
            print(msg)
            #raise IOError, msg
            print('upf will keep proceding with using upf.parse_epf method with n_unit=79')
            print('-'*52)
            blocks = parse_epf(filename)

        npf = len(blocks)
        if npf==0: raise IOError('No pf block found.')

        datasets = []; max_khi = []
        if  npf>1: hkls=[] ## multiple number of pole figures in a file
        for i in range(npf):
            #d = __epffiletrimmer__(blocks[i]) #only for popLA epf format
            ## epffiletrimmer should return a block of intensity grid map
            ## as well as maximum khi.
            ## --> modification (2011 NOV 8)
            hkl = blocks[i][0:3] #hkl
            hkl = list(map(int, [hkl[0],hkl[1],hkl[2]]))

            if npf>1: hkls.append(hkl)

            if blocks[i][3]!=')':
                print('Caution: unexpected hkl labeling format')

            d, mxk = __epffiletrimmer__(blocks[i])
            #only for popLA epf format
            datasets.append(d)
            max_khi.append(mxk)

        print("number of pole figures:", len(datasets))

        ## presumption of using 19x72 should be deprecated...
        data = np.zeros((len(datasets), 19, 72)) #npf, nphi, nkhi
        dum = data.copy()
        for i in range(len(datasets)):
            for j in range(len(datasets[i])):
                for k in range(len(datasets[i][j])):
                    data[i,j,k] = datasets[i][j][k]
        ## Swap the axis
        data = data.swapaxes(1,2) # from 72 x 19 to 19 x 72

        ## psuedo-normalization
        for i in range(len(data)):
            data[i] = pfnorm(data[i].copy())
        ## --
        if npf>1: hkl=hkls[::]
        else: hkl=[hkl]
        return data, max_khi, hkl

    elif mode=='xpc':
        """
        Adapted .xpc format parser, from
        https://github.com/usnistgov/texture
        commit 9c0ac85
        based on ready made popLA epf format parser
        """
        import pandas as pd
        print('You are now reading experimental pole figure(s) :%s'%filename)
        blocks = open(filename, 'rU').read().split('\n\n\n\n')[1:]
        print('There are %s blocks of data found'%len(blocks))
        if len(blocks)==0:
            msg1 = 'xpc parser in upf assumes that pole figures are separated by 4 new lines'
            msg2 = ' searching %s finds no set of 4 new lines in '%filename
            msg  = '%s \n %s'%(msg1,msg2)
            raise IOError(msg)
            # blocks = parse_epf(filename)

        npf = len(blocks)
        if npf==0: raise IOError('No pf block found.')

        datasets = []; max_khi = []
        if  npf>1: hkls=["HKL"] ## multiple number of pole figures in a file

        for part in blocks:
            line=part.split('\n')
            # print len(line)

            #header lines
            structureline=ff.FortranRecordReader('(6f10.4,1x,i4,1x,i4)') # the head-line in each data block
            [a,b,c,alpha,beta,gamma,crystalclass,something] \
                =structureline.read(line[1])
            pfDefline=ff.FortranRecordReader('(1x,3i3,6f5.1,2i2)')
            [h,k,l,unknown1,tilt,tiltinc,unknown2,rotation,\
             rotationinc,unknown3,unknown4]=pfDefline.read(line[2])

            #for the rest of the lines, do the following
            dataline=ff.FortranRecordReader('(1x,18i4)')

            # Pretty ugly code, but works...
            grouping=[[3,4,5,6],[7,8,9,10],[11,12,13,14],[15,16,17,18],[19,20,21,22],[23,24,25,26], \
                     [27,28,29,30],[31,32,33,34],[35,36,37,38],[39,40,41,42],[43,44,45,46],[47,48,49,50], \
                     [51,52,53,54],[55,56,57,58],[59,60,61,62],[63,64,65,66],[67,68,69,70],[71,72,73,74], \
                     [75,76,77,78]]

            #dataset=np.array([])
            dataset=[]
            for item in grouping:
                #print item[0],item[1],item[2],item[3]
                parsed=dataline.read(line[item[0]])
                parsed.extend(dataline.read(line[item[1]]))
                parsed.extend(dataline.read(line[item[2]]))
                parsed.extend(dataline.read(line[item[3]]))
                dataset.append(parsed)

            #print dataset

            #now saves as a dataframe, and wraps 360 to 0 degrees
            #row and column indexes are by degrees
            df=pd.DataFrame(dataset, index=np.arange(0,91,5))
            df.columns=[np.arange(0,360,5)]
            df[360]=df.ix[:,0]
            #print df.ix[:,0]

            #df["Tilt Angle"] =np.arange(0,92,5)
            #df.index.names=[np.transpose(np.arange(0,90,5))]
            #print df

            hkl = [h,k,l] #hkl
            #hkl = map(int, [hkl[0],hkl[1],hkl[2]])
            #print hkl
            hkls.append(hkl)
            datasets.append(df)
            #max_khi.append(mxk)

        #print hkls
        print("number of pole figures:", len(datasets))

        ### convert the panda data frame to numpy
        npf=len(datasets)
        datasets_compliant=np.zeros((npf,72,19))
        for i in range(npf):
            array = datasets[i].values
            arrayt = array.T
            ## ignore the last phi axis as it is repeated (phi360=phi0)
            datasets_compliant[i,:,:]=arrayt[:72,:]

        return datasets_compliant, hkls

    else: raise IOError('Unexpected mode is given')


def __epffiletrimmer__(block):
    """
    EPF (experimental pole figure) trimmer for a single block of data
    EPF: The popLA format

    Arguments
    =========
    block
    """
    pfdata = block.split('\n')
    dkhi = float(pfdata[0][5:9])   #khi incremental angle
    fkhi = float(pfdata[0][9:14])  #khi end angle
    dphi = float(pfdata[0][14:19]) #phi incremental angle

    # written as 360 but it is only upto 355
    fphi = float(pfdata[0][19:24]) - dphi

    ## now it is enfored to have 5 deg resolution, which
    ## results in 19x72
    ## intensities along khi and phi axes coordinate
    data = np.zeros((19,72))
    ## the rest information is neglected from the header

    pfdata = pfdata[1:]
    lines = []
    for i in range(len(pfdata)):
        if len(pfdata[i])==72:
            lines.append(pfdata[i][:])
        elif len(pfdata[i])==73:
            lines.append(pfdata[i][1:])

    pfdata = lines     #       (maxkhi/dkhi + 1) * 4
    if len(pfdata)!=76:# 76 =  ((90 / 5) + 1) * 4
                       # (4lines for one khi level)
        print('len(pfdata) =', len(pfdata))
        print(pfdata)
        raise IOError('Unexpected pfdata format or type')

    if True:
        for j in range(19): #90 / 5 + 1 #number of khi threads
            kh = pfdata[j*4: (j+1)*4] #along the same kh level
            khline = ''
            for k in range(4):
                khline = khline + kh[k]
            khline # all data in a thread of
                   # string every 4digits corresponds to a datum
            kp = 0 # each point along phi
            for k in range(72): #72 = 288 / 4 (4digits) : 0~355 5.0
                datum = khline[k*4:(k+1)*4]
                data[j,k] = int(datum)

    # block of mesh grid pole figure map and khi final angle.
    return data, fkhi

def shiftphi(data, dp, phi):
    """
    shifted phi modification

    Arguments
    =========
    data
    dp
    phi
    """
    tiny = 0.0000001
    if abs(phi[0] - 0.)> tiny and abs((dp-phi[0])-dp/2.)< tiny:
        dum = data.copy()
        for i in range(len(data)): #phi axis
            for j in range(len(data[i])): #khi axis
                dum[i, j] = (data[i-1, j] + data[i, j]) / 2.
        data = dum.copy()
    return phi, data

def duplicate0_360(phi, data):
    """
    If both phi=0 and phi=360 exist,
    delete phi=360

    Arguments
    =========
    phi
    data
    """
    tiny = 0.00000000001
    isallsame = True
    for i in range(len(data[0])):
        print(data[0,i])
        print(data[-1,i])
    if any(abs(data[0,i]-data[-1,i])>tiny for i in range(len(data[0]))):
        print(data[0,i])
        print(data[-1,i])
        isallsame = False

    isphioverlapped=False
    print(abs(phi[0] - phi[-1] - 2. * np.pi))
    if abs(phi[0] - phi[-1] + 2. * np.pi)<tiny:
        isphioverlapped=True

    if isphioverlapped and isallsame:
        newdata = np.zeros((data.shape[0]-1, data.shape[1]))
        for i in range(len(newdata)):
            for j in range(len(newdata[i])):
                newdata[i,j] = data[i,j]
        return phi[:-1], newdata
    elif isphioverlapped and isallsame==False:
        print("conflict results!")
        print("phi=0, and 360 turned out to coexist but")
        print("Values along these two axes are not equivalent")
        raise IOError
    elif not(isphioverlapped) and isallsame:
        print("Conflict results!")
        print("Phi is overlapped but phi[0] and phi[-1] is same")
        raise IOError
    else:
        print("No duplicated axis is found")
        return phi, data

def pole2f(poles_sa,poles_wgt,dth,dph,f):
    """
    Collect weights belonging to each 'sector' in the polar coordinate
    of pole figure projection.

    Arguments
    ---------
    poles_sa: pole referenced to sample axes
    poles_wgt: weight of each pole (originated from the weight of discrete orientations)
    dth: incremental rotation angle ( -180 ~  180)
    dph: incremental 'tilting' angle (   0 ~  180)
    """
    tiny = 1e-9
    for i, pole_sa in enumerate(poles_sa):
        theta,phi = cart2sph(pole_sa)
        ix = int((theta*180./np.pi+180)/dth-tiny)
        iy = int(  (phi*180./np.pi)    /dph-tiny)
        f[ix,iy]=f[ix,iy]+poles_wgt[i]
    return f
def pole2f_cols(poles_sa,poles_wgt,poles_col,f,dth,dph,fcol):
    """
    Collect weights belonging to each 'sector' in the polar coordinate
    of pole figure projection.

    Arguments
    ---------
    poles_sa: pole referenced to sample axes
    poles_wgt: weight of each pole (originated from the weight of discrete orientations)
    poles_col: extra column quantities
    f
    dth: incremental rotation angle ( -180 ~  180)
    dph: incremental 'tilting' angle (   0 ~  180)

    f_ij^{a}: weighted value of a scalar quantity {a}, belonging to each sector (ij).
    """
    tiny = 1e-9
    ncol=poles_col.shape[-1]
    for i, pole_sa in enumerate(poles_sa):
        theta,phi = cart2sph(pole_sa)#poles_sa[i])
        ix = int((theta*180./np.pi+180)/dth-tiny)
        iy = int(  (phi*180./np.pi)    /dph-tiny)
        for icol in range(ncol):
            fcol[ix,iy,icol]=fcol[ix,iy,icol]+poles_col[i,icol]*poles_wgt[i]/f[ix,iy]
    return fcol

def cart2polar(x,y):
    """
    cartesian to polar coordinate

    Arguments
    =========
    x
    y
    """
    r = np.sqrt(x**2+y**2)
    theta = np.arctan2(y,x)
    return r, theta

def circle(center=[0,0], r=1.):
    """
    Draw a circle around the given center point with
    a radius of r(as given)
    The default settings are as below.

    Arugments
    =========
    center = [0,0]
    r = 1.
    """
    ang = np.linspace(0,2*np.pi,1000)
    #unit circle * radius
    x = np.cos(ang)*r
    y = np.sin(ang)*r
    #circle transloation
    x = x + center[0]
    y = y + center[0]
    return x,y

def calc_vref_and_rot(a,b,csym,cdim,cang,nang=100):
    from .sym import cv

    icsym=get_icsym(csym)

    aca=cv(miller=a,icsym=icsym,cdim=cdim,cang=cang)
    bca=cv(miller=b,icsym=icsym,cdim=cdim,cang=cang)

    vref=np.cross(aca,bca)

    ##
    thf=np.arccos(np.dot(aca,bca))
    ths=np.linspace(0,thf,nang)
    rots=np.zeros((nang,3,3))
    varc=np.zeros((nang,3))

    for i, th in enumerate(ths):
        rots[i,:,:]=vector_ang(vref,np.rad2deg(th))
    return aca, bca, thf, vref, rots

# def get_ipf_boundary(a,b,c,nres,csym,cdim,cang):
#     pairs=[[a,b],[b,c],[c,a]]
#     coords=np.zeros((2,(nres-1)*3+1))

#     for i, pair in enumerate(pairs[:3]):
#         aca,bca,thf,vref,rots=calc_vref_and_rot(*pair,csym,cdim,cang,nres)
#         varc=calc_arc(aca,rots)
#         xy=np.zeros((len(varc),2))
#         for j, point in enumerate(varc):
#             xy[j,:]=projection(point)
#         i0=i*(nres-1)
#         i1=i0+nres-1
#         coords[:,i0:i1]=xy[0:-1,:].T
#     coords[:,-1]=coords[:,0]
#     return coords

def calc_arc(aca,rots):
    """
    Arguments
    ---------
    aca crystal direction
    rots: rotation matrix
    """
    nang=rots.shape[0]
    v_arc=np.zeros((nang,3))
    for i,rot in enumerate(rots):
        v_arc[i,:]=np.dot(rot,aca)
    return v_arc



def basic_triangle():
    """
    provide the boundary of basic triangle
    """
    v100_ca = np.array([1,0,0])
    v110_ca = np.array([1,1,0]) / np.sqrt(2.)
    v111_ca = np.array([1,1,1]) / np.sqrt(3.)

    delt_thet_100_110 = np.pi / 2.0
    delt_thet_111_110 =\
        np.arccos(np.dot(v110_ca, v111_ca))
    delt_thet_100_111 = \
        np.arccos(np.dot(v100_ca, v111_ca))

    line1 # 100_110
    line2 # 110_111
    line3 # 111_100


def trace(thet_f, n, v1, v2):
    """
    Trance of vectors that runs on the plane made
    by the two given vectors

    ---------
    Arguments
    ---------
    thet_f (in Radian)
    n
    v1
    v2
    """
    v1 = v1 / np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    v2 = v2 / np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    v3 = np.cross(v1, v2)
    v3 = v3 / np.sqrt(v3[0]**2 + v3[1]**2 + v3[2]**2)

    dthet = np.linspace(0., thet_f, n)
    # around v3 rotate v1 by incremental angle upto thet_f
    # rotate v1 about v3 axis.

    trace = []
    for iang in range(n):
        th = dthet[iang]
        r  = vector_ang(v3, th)
        vnew = np.dot(r, v1)
        trace.append(vnew)
    return trace

def xytrace(thet_f, n, v1, v2):
    """
    Returned projected trace
    """
    poles = trace(thet_f, n, v1, v2)
    X, Y = [], []
    for i in range(len(poles)):
        x,y = projection(poles[i], agrain=[0,0,0,1])
        X.append(x)
        Y.append(y)
    return X, Y

def vector_ang(u, th):
    """
    transformation matrix that rotates a vector about an axis u
    by the angle, th.

    ---------
    Arguments
    ---------
    u
    th (in Radian)
    """
    pi = np.pi
    idx = np.zeros((3,3))
    r   = np.zeros((3,3))
    for i in range(3):
        idx[i,i] = 1.

    ct = np.cos(th)
    st = np.sin(th)
    cm = crossop(u)

    for i in range(3):
        for j in range(3):
            r[i,j] = idx[i,j] * ct + st * cm[i,j] +\
                (1 - ct) * u[i] * u[j]

    return r

def crossop(u):
    m = np.zeros((3,3))
    m[0,1] = -u[2]
    m[0,2] =  u[1]
    m[1,0] =  u[2]
    m[1,2] = -u[0]
    m[2,0] = -u[1]
    m[2,1] =  u[0]
    return m

def __isunique__(a, b):
    """
    Is a(3) in b(m, 3)

    Arguments
    =========
    a
    b
    """
    for i in range(len(b)):
        ## is a//b[i]? either + or -
        diff0 = abs(b[i][0] - a[0]) + \
                abs(b[i][1] - a[1]) + abs(b[i][2] - a[2])
        diff1 = abs(b[i][0] + a[0]) + abs(b[i][1] + a[1])\
             + abs(b[i][2] + a[2])
        if diff0 < 0.1**4 or diff1 < 0.1**4:
            return False
    return True

def __circle__(center=[0,0], r=1.):
    """
    Draw a circle around the given center point with a
    radius of r(as given)
    The default settings are as below.

    Arugments
    =========
    center = [0,0]
    r = 1.
    """
    ang = np.linspace(0,2*np.pi,1000)
    #unit circle * radius
    x = np.cos(ang)*r
    y = np.sin(ang)*r
    #circle transloation
    x = x + center[0]
    y = y + center[0]
    return x,y

def deco_pf(ax,proj,triangle,cnt=None,miller=[0,0,0],
            iopt=0,iskip_last=False,
            ix='1',iy='2',mode='line',ires=True,
            nArray=None,levels=None,xcoord=None,ycoord=None,ilev=0,**kwargs_ipf):
    """
    Decorate matplotlib.pyplot.axes used for plotting
    (inverse) pole figures

    iopt==1: skip level lines

    Arguments
    ---------
    ax     (matplotlib axes)
    proj   ('pf' or 'ipf')
    triangle (* triangle boundary for inverse pole figure)
    cnt    (contour object generated by ax.contour)
    miller (miller indices)
    iopt
    iskip_last (whether or not to skip the maximum
                iso-contour line or level)
    ix     (xlabel, i.e., horizontal axis label)
    iy     (ylabel, i.e., vertial axis label)
    mode   pole figure plotting mode (line or fill)
    ires (logial) if ires true and mode is 'line'
         draw black dots on the background.
    nArray
    levels
    xcoord
    ycoord
    ilev (0) 0: 0: do not draw level lines 1: draw level lines
    **kwargs_ipf
    """
    from .sym import calc_cvec

    ## ------------------------------------------
    ## in case if inverse pole figures, obtain
    ## the size of triangle.
    if proj=='ipf':
        y1=triangle[1].max()
        y0=triangle[1].min()
        x1=triangle[1].max()
        x0=triangle[1].min()
        yscale=y1-y0
        xscale=x1-x0


    #--------------------------------------------
    ## fontsize of appended text to pole figures will be 4*fact
    fact = 1.8
    # --
    if mode in ['fill','line']:
        clev    = cnt._levels
        tcolors = cnt.tcolors
        if iskip_last: nlev = len(tcolors)-1
        else:  nlev = len(tcolors)

    #--------------------------------------------
    ## place colorbar of contours with its intensities
    if iopt==1: pass
    elif iopt==0 and mode in ['fill','line']:
        if proj=='pf':
            x=[1.35,1.39]
        elif proj=='ipf':
            x0=triangle[0].max()+0.02
            x1=x0+0.07
            x=[x0,x1]

        for i in range(nlev):
            cc = tcolors[i][0][0:3]

            if proj=='pf':
                y=[1. - i * 0.25, 1. - i * 0.25]
            elif proj=='ipf':
                y0=triangle[1].max()+0.2
                scale=y0/8
                y=[y0-i*scale, y0-i*scale]

            if not(iskip_last) and i==nlev-1 and mode=='line':
                ax.plot((x[0]+x[1])/2.,(y[0]+y[1])/2.,
                        '+',mew=2.,color=cc)
            else:
                ax.plot(x,y,color=cc)
            ## level text
            if clev[i]<10: s='  %4.2f'%clev[i]
            else:          s='%5.2f'%clev[i]
            ax.text(x=x[0]+0.02,
                    y=y[0],
                    s=s,fontsize=4.0*fact,va='center')

    #--------------------------------------------
    ## Place the Miller indices at the three
    ## corners of the triangle.
    if proj=='ipf':
        a=kwargs_ipf['a']
        b=kwargs_ipf['b']
        c=kwargs_ipf['c']
        fnsx=kwargs_ipf['fnsx']
        if ires and mode=='line':
            mask_invpf_tri=kwargs_ipf['mask_invpf_tri']
        loc=(0,0)
        for j, mil in enumerate([a,b,c]):
            t=''
            for i, v in enumerate(mil):
                if v<0: tx=r'\bar{%i}'%-v
                else: tx='%i'%v
                t=f'{t}%s'%tx
            t=rf'$({t})$'
            if j==0:
                x=triangle[0].min()-yscale/6.
                y=triangle[1].min()-yscale/6.
            if j==1:
                x=triangle[0].max()+xscale/6.
                y=triangle[1].min()-yscale/6.
            if j==2:
                x,y=projection(-calc_cvec(miller=mil,fnsx=fnsx))
                y=triangle[1].max()+yscale/6.+0.05
            loc=(x,y)
            ax.text(*loc,t,va='center',ha='center')

        ax.set_xlim(min(triangle[0])-0.05,max(triangle[0])+0.05)
        ax.set_ylim(min(triangle[1])-0.05,max(triangle[1])+0.05)

    #--------------------------------------------
    ## Add small block dots on the background
    ## in case contouring is done with lines
    ## but not filled.
    if ires and mode=='line' and ilev==1:
        filt=nArray[:,:]<levels[0]
        if proj=='ipf':
            filt=np.logical_and(filt,np.logical_not(
                mask_invpf_tri))
        filt[0,1:]=False
        filt[1:,0]=False

        filt=np.array(filt,dtype='bool')
        xs=xcoord[filt]; ys=ycoord[filt]

        if len(xs)>0:
            alpha=0.02
            #if ismooth>2: alpha=0.02
            ax.plot(xs,ys,'k.',alpha=alpha,
                    markersize=2.0)


    ##---------------------------------------------
    ## below lines are applied regardless of 'mode'.
    ax.set_axis_off()
    ax.set_aspect('equal')

    ##---------------------------------------------
    ## Annotate Miller indexed pole, either crystal
    ## or sample.
    s='('
    for k in range(len(miller)):
        if miller[k]<0: h = r'\bar{%s}'%str(-1*miller[k])
        else: h = '%s'%str(miller[k])
        s='%s%s'%(s,h)
    s='%s)'%s
    s=r'$\mathbf{%s}$'%s
    if proj=='pf':
        ax.text(0,-1.3,s,fontsize=9,ha='center',
                va='center')
    if proj=='ipf':
        scale=triangle[1].max()-triangle[1].min()
        ax.text(
            (min(triangle[0])+max(triangle[0]))/2.,
            min(triangle[1])-scale/3.,s,fontsize=9,
            ha='center',va='center')

    #----------------------------------------------#
    ## Triangle, circle & vertical and horizontal
    ## labels (ix, iy), and ticks.
    if proj=='pf':
        _x_,_y_ = __circle__()
        ax.plot(_x_,_y_,'k-')
        ax.set_xlim(-1.1,1.4)
        ax.set_ylim(-1.1,1.4)

        ## axis label/    ## Ticks
        ax.text(1.15,0. ,ix,va='center',ha='center')
        ax.text(0. ,1.15,iy,va='center',ha='center')
        ax.plot([0.0,0.0], [0.97,1.00],'k-')
        ax.plot([0.97,1.00],[0.0,0.0],'k-')

    if proj=='ipf':
        ax.plot(*triangle,'-k',zorder=1e10)


def projection(pole=None):
    """
    Projects a pole (vector) to projection plane.
    (default is stereography projection)

    pole = [1,1,1] or [1,1,0] something like this.

    The stereographic projection uses vectors pointing at
    southern hemisphere.

    Argument
    --------
    pole = None

    Returns
    -------
    X, Y
    """
    #normalization of the miller indices
    pole = pole / np.sqrt(pole[0]**2 + pole[1]**2 + pole[2]**2)
    a,b,c = pole[0:3]
    ###  mid-plane projection (z=0)
    if abs(c-1)<1e-8:
        # pole[0] = 0; pole[1]=0; pole[2] = 1
        X=0; Y=0
    else:
        X = a/(c-1)

        Y = b/(c-1)
    return X,Y

def invproj(x=None,y=None):
    """
    Converts the given projected point (x,y) to a pole in 3D.

    Arguments
    ---------
    x = None
    y = None

    Returns
    -------
    np.array([X,Y,Z])
    """
    X = 2*x/(1+x**2+y**2)
    Y = 2*y/(1+x**2+y**2)
    Z = (-1+x**2+y**2)/(1+x**2+y**2)
    return np.array([X,Y,Z])

def vect2sphe(pole):
    """
    cartesian vector to spherical cooridnate system

    Argument
    ========
    pole

    Returns
    -------
    x,y
    """
    seca = math.sqrt(pole[0]**2 + pole[1]**2)
    if seca < 0.1**6 :
        x = 0.0
        y = cos(seca)
    else:
        x = math.atan2(pole[1], pole[0])
        y = pole[2]
    return x, y

def cart2sph(pole):
    """
    argument: pole
    Returns phi and theta in such order indeed.

    Argument
    ========
    pole

    Returns
    -------
    [phi,theta]
    """
    pole = np.array(pole)
    r = np.sqrt((pole**2).sum())
    theta = math.acos(pole[2]/r)
    phi = math.atan2(pole[1],pole[0])
    return np.array([phi, theta])

def agr2pol(agrain=None, miller=None, proj=None):
    """
    -- for pole figure projection (proj='pf')
    For the given grain, crystallographic direction
    is mapped onto the sample axes by doing tensor
    product of the pole vector with rotation matrix.

    -- for inverse pole figure projection (proj='ipf')
    For the given grain, a vector in sa transformed to
    one referred in the crystal axes.

    Arguments
    =========
    agrain = None
    miller = None
    proj   = None
    """
    if proj==None: print("argument proj should be given"); raise IOError
    elif proj!='pf' and proj!='ipf':
        print(" proj should be either 'pf' or 'ipf'")
        raise IOError
    if type(miller).__name__=='list': miller = np.array(miller)

    # a-matrix between the grain's coordinate and sampl coordinate
    phi1 = agrain[0]; phi = agrain[1]; phi2 = agrain[2] #VPSC convention
    # rotation matrix (sa->ca)
    amat = euler(ph=phi1, th=phi, tm=phi2, echo=None)

    #normalize miller
    norm = math.sqrt( miller[0]**2 + miller[1]**2 + miller[2]**2)
    miller = miller / norm

    if proj=='pf':
        # Returns the dot product of the
        # transposed A-matrix and miller vector
        "A^{T}_{ij} * V_{j}"
        return np.dot(amat.transpose(), miller) #ca to sa

    elif proj=='ipf':
        # #for inverse projection,
        # #the sample axis should be one of principal axes.(x or y or z)
        # if   miller[0]==1: p_ca=0
        # elif miller[1]==1: p_ca=1
        # elif miller[2]==1: p_ca=2
        # # else: print"something is wrong"; raise IOError
        # # map the pole in the sample axes into the crystal axes
        "A_{ij} * V_{j}"
        return np.dot(amat, miller) #sa to ca returns the pole in ca
    else:
        print("projection should be pf or ipf")
        raise IOError

def ipfline(center=[0,0],csym='cubic'):
    """
    The boundary line for inverse pole figure
    ---------
    Arguments
    ---------
    center = [0,0]
    csym   = 'cubic'
    """
    xc = []; yc = []
    if csym!='cubic': print("Only Cubic!"); raise IOError
    xc.append( center[0])
    yc.append( center[1])

    for i in np.linspace(0.,1/math.sqrt(3.)):
        yaux = i
        xaux = math.sqrt((1. - yaux**2)/2.)
        zaux = xaux
        t1 = math.sqrt(1. - zaux) / math.sqrt(xaux**2 + yaux**2)
        t2 = t1/math.sqrt(1. + zaux)
        ## equal area
        # xc.append(xaux*t1)
        # yc.append(yaux*t1)
        ## stereo
        xc.append(xaux*t2)
        yc.append(yaux*t2)

    xc.append(center[0])
    yc.append(center[1])
    return np.array([xc,yc])

"""
Sample symmetry application is performed over RVE calculation
Refer to the RVE class in cmb.py module.
"""
def gen_fig(nrows=3,ncols=3,colsize=2.5,rowsize=2.,**kwargs):
    from matplotlib.gridspec import GridSpec
    gs=GridSpec(nrows=nrows,ncols=ncols)
    fig=plt.figure(figsize=(colsize*ncols,rowsize*nrows),**kwargs)
    axes=np.empty((nrows,ncols),dtype='object')
    for i in range(nrows):
        for j in range(ncols):
            ax=fig.add_subplot(gs[i,j])
            axes[i,j]=ax
    return fig,axes

class polefigure:
    # decides if the given set is in the texture file form or array
    def __init__(self, grains=None, filename=None, fnsx=None, csym=None,
                 ngrain=100, cdim=None, cang=None,
                 ssym=False, epf=None,epf_mode=None):
        """
        cdim=[1.,1.,1.6235] ## AZ31

        ----------------
        class polefigure
        ----------------
        Makes or accepts grains and saves them into the global
        variable self.gr. As a local method, there is 'def core'
        in which miller indices of poles (pole), are returned for
        a grain (agrain). The crystallographically equivalent poles
        are calculated if symmetry operation is allowed (isym).

        ---------
        Arguments
        ---------
        grains = None
        filename = None
        fnsx = None - Single crystal file used in VPSC-family codes
        csym = 'cubic' or'hexag'  #crystal symmetry
        ngrain = 100
        cdim=[1.,1.,1.]
        cang=[90.,90.,90.]
        ssym=False : Sample symmetry: not even started...
        epf = None : experimental pole figure file
        epf_mode = 'epf' or 'xpc'
        """
        from . import sym
        # The grain aggregte can be given either through a file or #
        # passing an array of them to the class directly.          #
        # either grains or filename                                #
        # if none of them is given, a 500-grains file is generated #
        # and returns its grains to the global gr variable.        #

        if type(grains)==type(None) and type(filename)==type(None)\
           and type(epf)==type(None):
            print(" ****************************** ")
            print(" Since no argument is passed,   ")
            print(" 1000 random grains are created ")
            print(" ****************************** \n")
            from .cmb import random
            self.gr = random(phi1=360,phi2=360,phi=90,ngrain=1000,iplot=False)

        self.epf = epf # global

        if type(grains)!=type(None):
            self.gr = np.array(grains)
        elif type(filename)!=type(None):
            with open(filename,'r') as fo:
                lines_ori=fo.readlines()
                lines=lines_ori[4:]

            ncol=len(lines[0].split())

            try:
                ## attempt to find ngr from the 4th line (works with 'TEX_PHx.OUT' format)
                # print(f'lines_ori[3]:, {lines_ori[3]}')
                ngr=int(lines_ori[3].split()[1])
            except:
                ngr=len(lines_ori)-4

            self.gr=np.zeros((ngr,ncol))
            lines=lines[-1:-1-ngr:-1][::-1]
            for i, line in enumerate(lines):
                self.gr[i,:]=np.fromiter(map(float,line.split()),float)
        elif type(epf)!=type(None): # None is the default for epf
            """
            experimental pole figures..
             # available formats:
                 - UXD
                 - steglich
                 - bruker
                 - epf*
                 - xpc
            """
            if type(epf).__name__=='list': self.epf_fn = epf
            elif type(epf).__name__=='str': self.epf_fn = [epf]
            elif type(epf)==type(True):
                fn = [] # list of file names
                print('type the experimental pole figure file names')
                print("To finish input, press enter")
                while True:
                    dum = input(">>> ")
                    if len(dum)==0: break
                    fn.append(dum)
                self.epf_fn = fn
            else: raise IOError('Unexpected epf type found')

            ## check if the file name is correct ##
            for i in range(len(self.epf_fn)):
                if not(os.path.isfile(self.epf_fn[i])):
                    raise IOError("Could not find %s"%self.epf_fn[i])
            ## --------------------------------- ##

            ## POLE FIGURE MODE --------------------------------------
            if type(epf_mode)==type(None):
                print("Type the experimental polfe figure mode")
                print("Available options:", end=' ') #continuation
                print("bruker, steglich, epf, xpc (default: %s)"%'epf')
                epf_mode = input(" >>>" )
                if len(epf_mode)==0: epf_mode='epf' # default
            ##---------------------------------------------------------

            self.grid = []; self.hkl = []
            ## more than one pf can be included.
            npole_per_file = []
            if epf_mode=='epf': self.max_khi = [] #Available only for epf_mode yet.

            for i in range(len(self.epf_fn)):
                if epf_mode=='epf':
                    data, maxk, hkl = epfformat(
                        mode=epf_mode,
                        filename=self.epf_fn[i])
                    # one file may include multiple poles
                    for i in range(len(data)):
                        self.grid.append(data[i])
                        self.max_khi.append(maxk[i])
                        self.hkl.append(hkl[i])
                    npole_per_file.append(len(data)) # of pole per a file
                elif epf_mode=='xpc':
                    data, hkl = epfformat(
                        mode=epf_mode,
                        filename=self.epf_fn[i])
                    for i in range(len(data)):
                        self.grid.append(data[i])
                        # self.max_khi.append(90.)
                        self.hkl.append(hkl[i])
                else:
                    data = epfformat(
                        mode=epf_mode,
                        filename=self.epf_fn[i])
                    self.grid.append(data)

                    self.hkl.append(None)
            self.grid = np.array(self.grid)
            self.epf_mode=epf_mode

        ## EXPERIMENTAL POLE FIGURE
        ## ------------------------------------------------------- ##
        ## POLE FIGURES BINNED FROM THE POLYCRYSTALLINE AGGREGATES ##

        if epf==None:
            dat = self.gr.transpose()
            phi1 = dat[0]; phi = dat[1]; phi2 = dat[2]
            ph1min, ph1max= int(
                round(min(dat[0]/90.)))*90, int(
                round(max(dat[0]/90.)))*90
            phmin, phmax  = int(
                round(min(dat[1]/90.)))*90, int(
                round(max(dat[1]/90.)))*90
            ph2min, ph2max= int(
                round(min(dat[2]/90.)))*90, int(
                round(max(dat[2]/90.)))*90

            ## symmetric multiplication over self.gr is performed unless ph1max==360
            """
            Sample symmetry application is pending,
            because it is done over rve (refer to cmb.py)
            """

            ### environments global variables
            #1 symmetry
            self.ngr = len(self.gr)

            if type(fnsx)!=type(None) and \
               (type(csym)!=type(None) or  \
                type(cdim)!=type(None) or \
                type(cang)!=type(None)):
                print('**Error: specify either fnsx or (csym,cdim,cang)')
                print('        But, do not specify both')
                print(f'fnsx: {fnsx}')
                print(f'csym: {csym}')
                print(f'cdim: {cdim}')
                print(f'cang: {cang}')
                raise IOError('** Error in fnsx/csym,cdim,cang')
            elif type(fnsx)==type(None) and \
               (type(csym)==type(None) or  \
                type(cdim)==type(None) or \
                type(cang)==type(None)):
                print('**Error: At least either fnsx or (csym,cdim,cang) should be given')
                raise IOError('** Error in fnsx/csym,cdim,cang')


            if type(fnsx)==type(None):
                self.fnsx = None
                self.csym = csym
                self.cdim = cdim
                self.cang = cang
            else:
                self.fnsx = fnsx
                self.csym, self.cdim, self.cang \
                    = sym.read_fnsx(self.fnsx)

    def epfplot(self,ifig,cmap,nlev,mn,mx,ix,iy,rot,iline_khi80):
        """
        This function is expected to be called
        within self.pf or self.pf_new

        if type(self.epf).__name__!='NoneType':
            ## This function is called within self.pf or self.pf_new

        * pole figure master data: self.grid
          [i]: i-th pole
          [i][j]: i-th pole's j-th phi segment's khi
          [i][j][k]: An intensity (count) at i-th pole's j-th phi sement at the k-th khi

        Arguments
        =========
        ifig
        cmap    ## color map, e.g., 'jet' or 'gray_r'
        nlev
        mn
        mx
        ix
        iy
        rot
        iline_khi80
        """
        from matplotlib.colors import LogNorm
        print('List of files:')
        for f in self.epf_fn: print('%s '%f)
        print('dimension of self.grid:', self.grid.shape)
        print('ifig=',ifig)

        fact = 2.
        nrow = len(self.grid)
        fig = plt.figure(figsize=(3.3*nrow,3.0))


        print('self.grid.shape:',self.grid.shape)
        mns,mxs,indices_mx = self.calcMXN(self.grid,mx,mn,'line',1)

        # print 'mns:',mns
        # print 'mxs:',mxs

        ## loop over each of pole figures
        for ip in range(len(self.grid)): #upon each of self.eps_fn
            ax = fig.add_subplot(1,nrow,ip+1)
            pf = np.zeros((self.grid[ip].shape[0]+1,
                           self.grid[ip].shape[1]))
            for i in range(len(self.grid[ip])):
                for j in range(len(self.grid[ip][i])):
                    pf[i,j] = self.grid[ip][i][j]
                    pf[-1,j] = self.grid[ip][0][j]

            # if type(mx).__name__=='NoneType':
            #     mx = np.array(pf).flatten().max()
            # if type(mn).__name__=='NoneType':
            #     mn = np.array(pf).flatten().min()

            # if mn==0: mn=0.5
            # if mx>100: mx=99.
            # levels = np.logspace(
            #     np.log10(mn),np.log10(mx),nlev)

            if mns[ip]==0: mns[ip]=0.5
            levels = np.logspace(
                np.log10(mns[ip]),np.log10(mxs[ip]),nlev)

            norm = LogNorm()

            nm = len(pf); nn = len(pf[0]) #phi, khi
            dp = 360. / nm; dk = 90. / nn
            khi = np.linspace(np.pi, np.pi/2., nn)
            phi = np.linspace(0., 2.*np.pi, nm)
            r   = np.sin(khi)/(1-np.cos(khi))
            R, PHI = np.meshgrid(r,phi)
            PHI = PHI + rot*np.pi/180. # default=0.

            x = R*np.cos(PHI); y = R*np.sin(PHI)

            print('levels in epfplot:')
            print(levels)

            pf[pf[::]<=0]=1e-9

            #cnt=ax.contour(
            cnt=ax.contourf(
                x,y,pf,levels=levels,
                cmap=cmap,norm=norm)
            deco_pf(ax=ax,cnt=cnt,miller=self.hkl[ip],ix=ix,iy=iy)

            ## dots like in pf_new.
            xs=[];ys=[]
            for j in range(len(x)-1):
                for k in range(len(x[j])):
                    if pf[j,k]<levels[0]:
                        if k==0 and j>1: pass
                        else:
                            xs.append(x[j][k])
                            ys.append(y[j][k])
            if len(xs)>0:
                ax.plot(xs,ys,'k.',
                        alpha=0.17,markersize=2.0)

            if iline_khi80:
                max_khi = 80.
                r_khi=2-np.sin(max_khi*np.pi/180.)/(1-np.cos(max_khi*np.pi/180.))
                rx,ry = __circle__(center=[0,0], r=r_khi)
                ax.plot(rx,ry,'--',color='gray')

        return fig

    def transformation(self,transf):
        """
        Apply transf to <self.gr>.

        Arguments
        ---------
        <transform>
           transformation matrix applied to the entire polycrystal aggregate.
        """
        for i in range(len(self.gr)):
            phi1,phi,phi2,wgt = self.gr[i][:4]
            ## amat = arg[-1]
            amat=euler(phi1,phi,phi2,a=None,echo=False) ## ca<-sa
            amat=amat.T ## sa<-ca
            if (transf==np.identity).all():
                pass
            else:
                amat=np.dot(transf,amat)

            phi1,phi2,phi3 = euler(a=amat.T)
            self.gr[i][:3]=phi1,phi2,phi3

    # def pf_axis(self, pole=[[1,0,0]], ifig=1):
    #     """
    #     Plot each pole without crystal symmetry
    #     """
    #     color = ['r','b','g','k','gray']
    #     #marker =['o','x','+','d','.']
    #     for ip in range(len(pole)):
    #         cl = color[ip]
    #         #mk = marker[i]
    #         for i in range(len(self.gr)):
    #             tm = self.dotplot(proj='pf', agrain=self.gr[i],
    #                               npole=len(pole), ipole=ip+1,
    #                               pole=pole[ip], ifig=ifig,
    #                               cdim='None', cang=self.cang,
    #                               csym=self.csym, mode=None,
    #                               color=cl)

    # def pf2xyw(self,pole=[1,0,0],csym='cubic',cdim=[1.,1.,1.],
    #            cang=[90.,90.,90.],fn='dat.xyz'):
    #     """
    #     Read pole and write xyw to a file

    #     Arguments
    #     =========
    #     pole
    #     """
    #     f = open(fn,'w')
    #     xyzw= []
    #     for i in range(len(self.gr)):
    #         gr = self.gr[i][::]
    #         phi1, phi, phi2 = gr[:3:]
    #         phi1 = phi1 - 90.

    #         npeq = __equiv__(
    #             miller=pole, csym=csym, cdim=cdim, cang=cang)

    #         xy, POLE = self.core(
    #             pole=pole, proj='pf',csym=csym,
    #             agrain=gr,isym=True,
    #             cdim=cdim,cang=cang, equivp=npeq)

    #         w = gr[-1]
    #         # for j in xrange(len(xy)):
    #         #     x,y = xy[j]
    #         #     z = 0
    #         #     f.write('%4.2f %4.2f %4.2f %11.4e\n'%(x,y,z,w))
    #         #     xyzw.append([x,y,z,w])


    #         for j in range(len(POLE)):
    #             xyz=POLE[j]
    #             x,y,z = xyz
    #             f.write('%4.2f %4.2f %4.2f %11.4e\n'%(x,y,z,w))
    #             xyzw.append([x,y,z,w])

    #     f.close()
    #     return np.array(xyzw).T

    def ipf(self, pole=None,ifig=4,mode='dot',deco=True,**kwargs):
        """
        Given the pole plot the dot inverse pole figure.
        **The contour version of the inverse pole
        figure is not available yet.

        ---------
        Arguments
        ---------
        pole  = [1,0,0]
        ifig  = 1
        mode  = 'dot', 'contour'
        deco  = True
        **kwargs - key-worded argument passed to plots.
        """
        if type(pole)==type(None):
            raise IOError('miller index of the pole should be given')

        if mode=='dot':
            temp = []
            for i in range(len(self.gr)):
                if deco==True and i==0:
                    _deco_=True
                else:
                    _deco_=False

                tm, fig = self.dotplot(proj='ipf',csym=self.csym,
                                       agrain=self.gr[i],
                                       pole=pole, ifig=ifig, deco=_deco_,
                                       **kwargs)
                temp.append(tm)
            #self.dotplot(proj='ipf',
            #agrain=self.gr[i], pole=[0,1,0], ifig=5)
            #self.dotplot(proj='ipf',agrain=self.gr[i],
            #pole=[0,0,1], ifig=6)
            return temp
        elif mode=='contour':
            fig    = plt.figure(ifig)
            ipf_ax = fig.add_subplot(111)
            for i in range(len(self.gr)):
                pass
            pass
        else: raise IOError('Unexpected model for ipf')

    def core(self, pole=None, proj='pf', csym=None, ssym=None,
             agrain=None, isym=True, cdim=[1,1,1],
             cang=[90.,90.,90.],
             equivp=None):
        """
        --------------------------------
        The core of the polefigure class
        --------------------------------

        It is the engine for plotting regular and inverse
        pole figures by generating projected
        cartesian coordinates for the 3D vectors
        onto the pole figure sphere (stereographic sphere).
        One can directly plot pole figure on a xy plane.

        Provided the miller indices,
        1. Calculates the crystallographically equivalent poles
        2. Maps the miller indices for a grain
           into a vector in sample axes
        3. Returns mapped poles as raw and as projected (xy)

        ---------
        Arguments
        ---------
          pole = None
          proj = 'pf'
          csym = 'cubic', 'hexag'
          ssym = None,
          agrain = None
          isym = True
          cdim = [ 1., 1., 1.]
          cang = [90.,90.,90.]
          equivp = None  #crystallographically equivalent pole
        """
        xy = []; POLE = [];
        if csym!='cubic' and csym!='hexag' and csym!='None' \
                and csym!='centro':
            raise IOError("Other symmetries than cubic"+\
                " or hexag or 'None'"+" nor 'centro' is"+\
                " not prepared yet")
        if proj!='pf' and proj!='ipf':
            raise IOError("Other modes of projection than pf and"+\
                "ipf is not prepared yet")

        if type(agrain)==type(None):
            raise IOError("A grains must be given to the method")
        if type(pole)==type(None):
            raise IOError("Pole must be given to core")

        if type(pole).__name__=='ndarray': pass
        elif type(pole).__name__=='list': pole = np.array(pole)

        else: raise IOError('Unexpected type of the pole argument')

        temp = pole.copy()
        del pole; pole = temp

        ##### calculates crystallographically equivalent poles #####
        ## pole figure projection
        if proj=='pf':
            if isym: npeq = equivp # xtallographically equiv. poles
            else: npeq = [pole]
            nit = len(npeq)
            nppp = []
            for i in range(nit):
                nppp.append(npeq[i])
                nppp.append(npeq[i]*-1)
            npeq = np.array(nppp)
            for i in range(len(npeq)):
                for j in range(len(npeq[i])):
                    if abs(npeq[i,j])<1e-9:
                        npeq[i,j] = 0.

        ## inverse pole figure
        elif proj=='ipf':
            # if abs(pole[0]**2+pole[1]**2+pole[2]**2-1)>0.1**6:
            #     print "The pole must be one of principal axes of"\
            #         " the sample"
            #     print "It should be among [1,0,0], [0,1,0], [0,0,1]"
            #     print "current pole is as below\n", pole
            #     raw_input()
            #     raise IOError
            npeq = [pole] ## Note it is sample axis vector!
            #equivalent pole calculatation is deffered to
            #next for in block
        ## unexpected proj argument
        else: print("It should be either pf or ipf"); raise IOError

        t0=time.time()
        t_agr2pol = 0.
        t_proj = 0.


        for ip in range(len(npeq)):
            ## 'pf':  converts ca pole to sa pole
            ## 'ipf': converts sa pole to ca pole
            t_1 = time.time()
            if   proj=='ipf':
                p = agr2pol(agrain=agrain, miller=npeq[ip], proj=proj)
            elif proj=='pf':
                ag = agrain[:3].copy()
                _np_eq_=npeq[ip]
                p = agr2pol(ag, _np_eq_,proj=proj)
            # p = agr2pol(agrain=agrain, miller=npeq[ip], proj=proj)

            t_agr2pol = t_agr2pol + time.time()-t_1
            if proj=='pf': # p is in sa
                ## if a pole is toward the north pole of the unit circle,
                #if p[2]>0: pass
                #else:
                POLE.append(p)
                t_1=time.time()


                xy.append(projection(pole=p))

                t_proj = t_proj + time.time()-t_1
            elif proj=='ipf': # p is in ca
                ## calculates equivalent p by
                ## applying symmetry operations
                if isym:
                    """
                    Sample axis is now referred to crystal
                    coordinate system. That vector has to
                    be mutiplicated by symmetry operations.
                    """
                    npoles = __equiv__(
                        miller=p, csym=csym,
                        cdim=cdim, cang=cang)
                    temp = []
                    for i in range(len(npoles)):
                        temp.append(npoles[i])
                        temp.append(npoles[i]*-1)
                    temp = np.array(temp)
                    npoles = temp

                else: npoles=[p]
                for npp in range(len(npoles)):

                    prj_xy = projection(pole=npoles[npp])
                    xy.append(prj_xy)
                    POLE.append(npoles[npp])
                pass # if over 'pf' or 'ipf'
            pass # End of for loop over ipf

        # print 'Elapsed time for agr2pol ', t2s(t_agr2pol)
        # print 'Elapsed time for proj ', t2s(t_proj)

        return xy, POLE

    def pf_new(self,ifig=None,axs=None,proj='pf',poles=[[1,0,0],[1,1,0]],ix='1',iy='2',
               mode='line',dth=10,dph=10,n_rim=2,cdim=None,ires=True,
               mn=None,mx=None,lev_norm_log=True,nlev=7,ilev=1,levels=None,
               cmap='magma',rot=0.,iline_khi80=False,
               transform=np.array([[-1,0,0],[0,-1,0],[0,0,1]]),
               ideco_lev=True,ismooth=1,
               **kwargs):
        """
        New version of pf that will succeed upf.polefigure.pf
        Note that upf.polefigure.pf is deprecated and will be deleted soon.

        Arguments
        ---------
        <ifig> or <axs>
            <ifig> and <axs> should be mutually exclusive.
            It is acceptable for both to be *not* specified.
            However, it is unacceptable for both to be specified.
        <proj>
           proj can be either 'pf' or 'ipf'. The former is the usual
           pole figure, while 'ipf' refers to the inverse pole figure
        <poles>
           For cubics, three digits; for hexagonals four digits
        <ix>, <iy>
           x and y tick labels appended to each pole figure
        <dph>:
            Grid of tilting angle
        <dth>:
            Grid of in-plane rotation angle
        <rot>:
             in-plane rotatation (radian)
        <n_rim>:
             The number of 'central' rims to be *averaged*.
             For better understandings, see the algorithm notebook
             located in ./ipynb/UPF_Algorithm.ipynb
        <cdim>:  crystal dimension
           For cubic, [1,1,1]
           For a particularly AZ31 sheet, it is [1,1,1.6235]
           Users should know what is the lattice dimension for the crystal
           structure of which he/she plots the pole figures.
        <ires>  = True;
           If True, indicate the grid
           If <mode> is 'fill' and ires is True, overlay the resolution
              all over the pole.
           if <mode> is 'line' and ires is True, only the spots lower than
              minimum level of contour is plotted.
        <mn>,<mx>
          minimum and maximum levels of contour
          If not specified, mn and mx is determined using levels of
          calculated contours
          if <mode> is 'fill', <mn> is overriden by the levels of
          calculated contours.
        <lev_norm_log>
           If True, use logarithmic scales. If False, linear scale.
        <nlev> = 7
           Level of iso contour bins.
           The total number of lines will be nlev+1
        <cmap>
           Color map used to color-code the contour levels.
           Refer to 'http://matplotlib.org/users/colormaps.html'
           for more color-map options.
        <iline_khi80> = False
           Whether or not to draw a line of chi=80 in experimental
           pole figure plot. - I usually obtain incomplete pole figure
           upto a tilting <chi> of 80.
        <mode>
           Contour modes: 'line', 'contour', or 'fill'
           dot modes    : 'dot' or 'dotm', or 'dotc'
             ** The option 'dotm' provides the coordinates and quits)
        <ilev>
           level option: 0 common contour levels for all poles generated
                         1 individual levels applied for individual poles
        <levels>
           Default is None. One can define levels of the contours.
        <transform>
           transformation matrix applied to the entire polycrystal aggregate.
        <ideco_lev> True or False
           switch to turn on or off the levels
        <ismooth>=1

        Returns
        -------
        fig: matplotlib.figure.Figure
        """
        import scipy

        ## check mutually exclusive arguments (ifig and axs)
        if type(ifig)!=type(None) and type(axs)!=type(None):
            raise IOError('** Err: ifig and axs are mutually exclusive')

        ####################################################
        ## PF plots for experimental pole figure is
        ## separately conducted by epfplot function
        if type(self.epf).__name__!='NoneType':
            print('** Writing Experimental pole figures')
            if transform!=np.identity(3).all():
                print('<transform> is ignored when plotting EPF')
            return self.epfplot(
                ifig=ifig,cmap=cmap,nlev=nlev, mn=mn, mx=mx,
                ix=ix,iy=iy,rot=rot,iline_khi80=iline_khi80)

        nlev = nlev + 1 ##
        miller=poles[::]
        if type(cdim)!=type(None): self.cdim=cdim

        if self.gr.shape[-1]>4: Ncol=[] ## esgr format

        ####################################################
        ## contoured or dotted (inverse) pole figures.

        #----------------------------------------------------
        ## Either create new matplotlib figure object or use
        #  the given <ifig> or even specifical  matplotlib
        #  axes objects.
        if type(axs)==type(None):
            if type(ifig)==type(None):
                fig = plt.figure(figsize=(3.3*len(poles),3.0))
            else:
                fig = plt.figure(ifig,figsize=(3.3*len(poles),3.0))
            ##
            axs=np.empty(len(poles),dtype='object')
            for i in range(len(poles)):
                axs[i] = fig.add_subplot(1,len(poles),i+1)
            plt.subplots_adjust(left=0,right=0.8)

        #--------------------------------------------
        ## for inverse pole figure, obtain the fund.
        ## triangle, and the three corner poles
        ## i.e., a,b, and c, as well as the masking
        ## array to remove contours outside the
        ## triangle region.
        triangle=None
        if proj=='ipf':
            ## stereographic triangle boundary
            triangle,a,b,c=get_ipf_boundary(fnsx=self.fnsx,nres=30)

        #----------------------------------------------------
        # contour (line, contour, fill) or dots (dot, dotm)
        if mode in ['line','contour','fill']:
            #------------------------------------------------
            # Calculate pole figure intensity nodes for
            # contouring.
            t0=time.time()
            N=[]
            for i in range(len(poles)):
                rst=cells_pf(
                    0,proj=proj,pole=poles[i],dth=dth,dph=dph,
                    csym=self.csym,cang=self.cang,
                    cdim=self.cdim,grains=self.gr,
                    n_rim = n_rim,transform=transform)
                if self.gr.shape[-1]>4: ## if extra columns exist
                    N.append(rst[0])
                    Ncol.append(rst[1])
                else:
                    N.append(rst)
            nArray=np.array(N)

            et = time.time()-t0
            try: uet(et,head='Elapsed time for calling cells_pf')
            except:pass

            #------------------------------------------------
            ## Get meshgrids
            R,Phi,x,y=get_grid_from_angles(dth,dph,rot)
            x_ori=x[::]; y_ori=y[::]

            #------------------------------------------------
            ## Smoothing by using scipy.ndimage.zoom
            if ismooth>1:
                t0=time.time()
                N_smooth=[]
                for i in range(len(poles)):
                    refined=scipy.ndimage.zoom(
                        nArray[i].T,ismooth)
                    # Remove unwanted negative intensities
                    refined[refined<0]=0.
                    N_smooth.append(refined.T)
                    # Newly refined (x,y) coordinates
                    # (required only once)
                    if i==0:
                        x,y=xy_grid_from_shape(refined.T.shape)
                nArray=np.array(N_smooth)

            #------------------------------------------------
            # Obtain the maximum and minimum intensities
            # Also, find the location of 'maximum' pole by
            # its index.
            mns, mxs, indices_mx = self.calcMXN(
                nArray,mx,mn,mode,ilev)

        #----------------------------------------------------
        # dotted pole figures
        elif mode in ['dot','dotm','dotc']:
            pf_dots=[]
            pf_dots_wgt=[]
            pf_dots_col=[]
            for i in range(len(poles)): # crystal or sample poles
                print(f'poles[i]:: {poles[i]}')
                XY,wgt,col_val=cells_pf(
                    1,proj=proj,pole=poles[i],dth=dth,dph=dph,
                    csym=self.csym,cang=self.cang,
                    cdim=self.cdim,grains=self.gr,
                    n_rim = n_rim,transform=transform)

                ## masking XY ...
                if proj=='ipf':
                    XY,tags=get_within_triangle(triangle,XY)
                    # use this tags to trim.
                    wgt=wgt[tags]
                    col_val=col_val[tags,:]

                pf_dots.append(XY)
                pf_dots_wgt.append(wgt)
                pf_dots_col.append(col_val)

            pf_dots=np.array(pf_dots,dtype='object')
            pf_dots_wgt=np.array(pf_dots_wgt,dtype='object')
            pf_dots_col=np.array(pf_dots_col,dtype='object')

            if mode=='dotm':
                return pf_dots
        else:
            raise IOError('Unexpected mode given to pf_new')

        #--------------------------------------------
        ## Create masking array to remove contours outside the
        ## triangle region.
        triangle=None
        mask_invpf_tri=None
        if proj=='ipf':
            ## stereographic triangle boundary
            triangle,a,b,c=get_ipf_boundary(fnsx=self.fnsx,nres=30)
            # if mode in ['line','contour','fill']:
            ## Generate masks to hide contours outside of the triangle
            if mode in ['line','contour','fill']:
                mask_invpf_tri=gen_mask_contour(
                    triangle,
                    shape=x.shape,x=x,y=y)

        #----------------------------------------------------
        # decoarting (inverse) pole figures
        # levels, boundaries (circle or triangle)
        # coloring of contours, miller indices for each poles
        for i in range(len(poles)):

            #------------------------------------------------
            ## determine the plotting function
            if   mode=='line':
                func = axs[i].contour
            elif mode in ['fill', 'contour']:
                func = axs[i].contourf
            elif mode in ['dot','dotc']:
                func = axs[i].scatter

            #------------------------------------------------
            # Decorating the contoured (inverse) pole figures
            if mode in ['line','contour','fill']:

                #--------------------------------------------
                ## create necessary objects for color-contours
                levels, cm_norm, cmap_mpl, color_mapping =\
                    get_pf_color_map_and_norm(
                        levels,cmap,lev_norm_log, mns[i],
                        mxs[i], nlev)

                #--------------------------------------------
                ## contour plot
                nArray[i][np.isnan(nArray[i])]=0. ## remove nan
                nArray[i][nArray[i]<=0]=1e-4      ## remove negative values.
                if proj=='ipf':
                    rst_within_triangle=np.ma.array(nArray[i],mask=mask_invpf_tri)
                    cnts=func(x,y,rst_within_triangle,levels=levels,
                              cmap=cmap,norm=cm_norm,zorder=10)
                else:
                    cnts=func(x,y,nArray[i],levels=levels,
                              cmap=cmap,norm=cm_norm,zorder=10)

                #--------------------------------------------
                # Indicate maximum location only for 'line' contours
                # but not the 'filled' contours.
                if mode=='line':
                    ## x, y coordinates of maximum intensity in grid
                    i0,j0 = indices_mx[i]
                    mx_coord_x = x[i0,j0]
                    mx_coord_y = y[i0,j0]
                    axs[i].plot(mx_coord_x,mx_coord_y,'+',mew=2,
                                zorder=1e2,
                                color=color_mapping.to_rgba(levels[-1]))

                #--------------------------------------------
                ## decorating (inverse) pole figures
                # if ideco_lev:ideco_opt=0
                # elif ideco_lev==False and i==len(poles)-1: ideco_opt=1
                # else: ideco_opt=0
                if ideco_lev and i==len(poles)-1: ideco_opt=0
                else: ideco_opt=1

                #if (ilev==1 or (ilev==0 and (i==len(axs)-1) or i==len(axs)-2)):
                ## arguments commonly used for PF and IPF
                kws=dict(ax=axs[i],proj=proj,
                         triangle=triangle,cnt=cnts,
                         miller=miller[i],iopt=ideco_opt,
                         iskip_last=False,ix=ix,iy=iy,
                         mode=mode,ires=ires,nArray=nArray[i,:,:],
                         levels=levels,xcoord=x,ycoord=y,ilev=ilev)

                if proj=='ipf':
                    kws.update(a=a,b=b,c=c,fnsx=self.fnsx,
                               mask_invpf_tri=mask_invpf_tri)
                deco_pf(**kws)

            #------------------------------------------------
            # Decorating the dotted (inverse) pole figures
            elif mode in ['dot','dotc']:
                if ideco_lev: ideco_opt=0
                else: ideco_opt=1
                x=pf_dots[i][:,0]
                y=pf_dots[i][:,1]

                if mode=='dotc':
                    ## color
                    ## create necessary objects for color-contours
                    levels, cm_norm, cmap_mpl, color_mapping =\
                        get_pf_color_map_and_norm(
                            levels,cmap,lev_norm_log,
                            pf_dots_wgt[i].min(),
                            pf_dots_wgt[i].max(),nlev)

                    print('min and max wgt:', pf_dots_wgt[i].min(),
                          pf_dots_wgt[i].max())

                    kwargs.update(c=pf_dots_wgt[i])

                func(x,y,**kwargs)

                kws=dict(ax=axs[i],proj=proj,
                         triangle=triangle,cnt=None,
                         miller=miller[i],iopt=ideco_opt,
                         iskip_last=False,ix=ix,iy=iy,
                         mode=mode,ires=ires,nArray=None,
                         levels=None,xcoord=x,ycoord=y)

                if proj=='ipf':
                    kws.update(a=a,b=b,c=c,fnsx=self.fnsx,
                               mask_invpf_tri=mask_invpf_tri)
                deco_pf(**kws)

        #--------------------------------------------------#
        ## returning some objects.
        if self.gr.shape[-1]>4 and proj=='pf':
            return fig, np.array(N), np.array(Ncol), \
                R*np.cos(PHI),  R*np.sin(PHI)
        elif self.gr.shape[-1]==4:
            try: return fig
            except: pass
        #--------------------------------------------------#





    def calcMXN(self,nArray=None,mx=None,mn=None,mode='line',ilev=0):
        """
        Calculate minimum and maximum

        Arguments
        ---------
        nArray: ndarray that contains all pole figure nodes.
        mx
        mn
        mode
        ilev  : option (0: common mx and mn, 1: individual mx and mn)

        Returns
        -------
        mns, mxs, indices_mx
        """
        npole = nArray.shape[0]
        mxs=np.zeros(npole)
        mns=np.zeros(npole)

        if ilev==0:
            ## determine maximum and minimum levels.
            if type(mx)==type(None): mx = nArray.flatten().max()
            if mx>100: mx=99.

            if type(mn)==type(None) and mode!='fill': mn = 0.5
            elif type(mn).__name__=='float':
                pass
            else:
                mn = nArray.flatten().min()
            ## commonly assigned
            mxs[:] = mx
            mns[:] = mn
        elif ilev==1:
            for ipole in range(npole):
                if type(mx)==type(None):
                    mx_ = nArray[ipole].flatten().max()
                else:
                    mx_ = mx*1.

                if mx_>100: mx_ = 99.

                if type(mn)==type(None) and mode!='fill':
                    mn_ = nArray[ipole].flatten().min()
                elif type(mn)==type(None) and mode=='fill':
                    mn_ = nArray[ipole].flatten().min()
                else:
                    mn_ = mn*1.

                mxs[ipole]=mx_
                mns[ipole]=mn_
        else:
            raise IOError

        ## Find coordinates of 'maximum' intensity
        indices_mx=[]
        for ipole in range(npole):
            indx_mx=np.unravel_index(nArray[ipole,:,:].argmax(),
                                     nArray[ipole,:,:].shape)
            ## indx_mx = np.argmax(nArray[ipole,:,:],axis=1)
            indices_mx.append(indx_mx)
        indices_mx=np.array(indices_mx)
        return mns, mxs, indices_mx

    def dotplot(self, pole=None, ifig=None, npole=1,
                ipole=1,
                proj='ipf', csym='cubic',
                ssym='tric', agrain=None,
                isym=True, irdc=True,
                cdim=[1.,1.,1.],
                cang=[90.,90.,90.], mode='trace', # or None
                deco=True,
                **kwargs
                # alpha=1.0, color='k',
                ):
        """
        Plots the individual poles. No contour.

        ---------
        Arguments
        ---------
        pole  = None
        ifig  = None
        npole = 1
        ipole = 1
        proj = 'pf' or 'ipf'
        csym = 'cubic' Crystal symmetry
        ssym = 'tric'  Sample symmetry
        agrain = None
        isym = True:   Symmetry operation to the pole
        irdc = True:   Reduction of inverse pole figure region
        cdim = [1., 1., 1.]
        cang = [90.,90.,90.]
        mode = 'trace'
        deco = True
        **kwargs : matplotlib pyplot key-worded arguments
        """
        if pole==None:
            print("Pole must be given")
            raise IOError

        ## If a grain is assigned,
        ## Operations are carried on it, instead of on
        ## Global grains.
        # if color==None:
        #     print 'Default color is black'; color='k'
        if type(agrain)!=type(None):
            "* a grains' euler angle is given"
            "* dot plotting is performed on a grain"
            agr = np.array([agrain]) #p1,p,p2, inten
        else:
            print('Agrain must be asigned')
            raise IOError
        #print "------------------------------------"
        XY = []; polz = []
        ## Draws the boundary of inverse and normal pole figures
        if type(ifig)==type(None):
            fig=None
            pass
        else:
            ## if plt.figure(ifig) is present, use it
            ## otherwise, create one.

            fact = 3.
            figsize = (npole*fact, 1.*fact)
            if plt.fignum_exists(ifig):
                fig=plt.figure(ifig)
            else:
                fig = plt.figure(ifig, figsize=figsize)

            ax  = fig.add_subplot(1,npole,ipole)

            if deco:
                ax.set_axis_off(); ax.set_aspect('equal')
                ax.text(x=-0.08, y=-0.07, s='(100)', fontsize=4.*fact, transform=ax.transAxes)
                ax.text(x=0.7, y= -0.07, s='(110)', fontsize=4.*fact, transform=ax.transAxes)
                ax.text(x=0.65, y= 0.8, s='(111)', fontsize=4.*fact, transform=ax.transAxes)
                # ax.text(x=0.5, y=-1.05, s='(%i%i%i)'%
                #         (pole[0],pole[1],pole[2]), fontsize=3. * fact,transform=ax.transAxes)
            if proj=='pf' and deco:
                rx, ry = __circle__(center=[0,0], r=1)
                ax.plot(rx, ry, color='grey')
                ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2)
            elif proj=='ipf' and deco:
                # in ipfline module, projection can be
                # changed into equal area type
                cxy = ipfline(center=[0,0], csym=csym)
                ax.plot(cxy[0], cxy[1], color='grey', alpha=0.1)

        ## add poles
        for i in range(len(agr)):
            ### crystallographically equivalent poles are calculated.
            agr[i][0] = agr[i][0] - 90.
            npeq = __equiv__(
                miller=pole, csym=csym, cdim=cdim, cang=cang)
            ### -----------------------------------------------------
            xy, POLE = self.core(
                pole=pole, proj=proj, csym=csym,
                agrain=agr[i], isym=isym, cdim=cdim,
                cang=cang,equivp = npeq)

            ##### POLE FIGURE PROJECTIN #####
            if proj=='pf':
                for j in range(len(xy)):
                    if POLE[j][2]<=0: #north pole is [0,0,1]
                        XY.append(
                            [xy[j][0], xy[j][1], agr[i][3]]
                            ) #x,y, intensity

            ##### INVERSE POLE FIGURE PROJECTION #####
            elif proj=='ipf':
                for j in range(len(xy)):
                    r = np.sqrt(xy[j][0]**2 + xy[j][1]**2)
                    if r>1.0: pass
                    else:
                        ### reduced region filter
                        ### must be applied here.
                        # y must be positive.
                        tiny = 0.
                        phi = math.atan2(xy[j][1],xy[j][0])
                        phi = phi * 180.0 / math.pi
                        if phi > 0. - tiny and phi < 45.0 +tiny:
                            ## another filter
                            ## 2nd and 3rd triangles
                            #a = POLE[j][0]; c = POLE[j][2]
                            a,b,c = invproj(x=xy[j][0], y=xy[j][1])
                            #print 'a,b,c= ',a,b,c
                            #print 'atan2(a,c)',
                            #math.atan2(a,-c)*180./math.pi
                            if math.atan2(
                                a,-c)*180./math.pi < 45.0 + tiny:
                                XY.append([xy[j][0], xy[j][1],
                                           agr[i][3]])
                                   #x, y, intensity of the grain

        if type(ifig)!=type(None):
            # plotting --------------------
            try:
                XY = np.array(XY)
                xxx = XY.copy().transpose()[0]
                yyy = XY.copy().transpose()[1]
                ax.plot(xxx,yyy,ls='None', **kwargs)
                if proj=='pf':
                    ax.set_xlim(-1.1,1.1)
                    ax.set_ylim(-1.1,1.1)
                # elif proj=='ipf':
                #     ax.set_xlim(-0.02,0.5)
                #     ax.set_ylim(-0.02,0.5)
            except:
                if len(XY)==0:
                    pass
                    #   print 'Empty XY is returned'
                else:
                    raise IOError("Unexpected Error");
                    # raw_input()
            #-------------------------------
        return np.array(XY),fig

def cells_pf(iopt=0,proj='pf',pole=[1,0,0],dph=7.5,dth=7.5,csym=None,cang=[90.,90.,90.],
             cdim=[1.,1.,1.],grains=None,n_rim=2,transform=np.identity(3)):
    """
    For the case of contours, creates cells gridded in
    the dimension of (nphi, ntheta). Given the delta x
    and delt y (dm, dn), each pole's weight is assigned
    to the cell to which it belongs.

    Plots the cell's weight and returns the cell in array.

    + Additional feature (2024-05-28)
    Some more updates that I'm trying to implement here is to visualize
    some other 'scalar' quantities that come along with the form of discrete orientation files.
    Conventionally, 'texture' file has only 4 columns, consisting of phi1, phi, phi2, and weight.
    The case that I'd like to deal here is to use the additional columns denoting, say, stored
    energy, in the form of pole figure.

    + Additional feature (2024-06)
    I plan to use this module for both pole figures and inverse pole figures

    ---------
    Arguments
    ---------
    <iopt> (0: contouring; 1: dots)
    <proj> can be either 'pf' or 'ipf'
    <pole> = [1,0,0]
    <dph>  = 7.5. (tilting angle : semi-sphere 0, +90 or full-sphere 0, +180)
    <dth>  = 7.5. (rotation angle: -180,+180)
    <csym> = None
    <cang> = [90.,90.,90.]
    <grains> = None, [] array shape: (ngr, 3)
    <n_rim>=2
    <transform>: the default is np.identity(3)

    ---------
    Returns
    -------
    nodes             if ncols==4  (just like the texture file of (E)VPSC)
    nodes, nodes_col  if ncols>4
    """
    from . import sym

    tiny = 1e-9
    ## Set of equivalent vectors based on crystal symmetry
    pole=np.array(pole)
    if proj=='pf':
        p0 = __equiv__(miller=pole,csym=csym,cdim=cdim,cang=cang)
        poles_ca=np.zeros((len(p0)*2,3))
        # both forward and backward poles
        poles_ca[:len(p0),:] = p0[:,:] # (+)
        poles_ca[len(p0):,:] =-p0[:,:] # (-)
        poles_ca         = poles_ca / np.sqrt((poles_ca**2).sum())
        ## poles_projected can be either crystal poles (PF) or sample poles (IPF)
        poles_projected  = np.zeros((len(grains),len(poles_ca),3))
        poles_wgt        = np.zeros((len(grains),len(poles_ca)))
        poles_col        = np.zeros((len(grains),len(poles_ca),grains.shape[-1]-4))

        for i, gr in enumerate(grains):
            phi1,phi,phi2,wgt = gr[:4]
            amat=euler(phi1,phi,phi2,a=None,echo=False) ## ca<-sa
            amat=amat.T ## sa<-ca
            if (transform==np.identity(3)).all():pass
            else:amat=np.dot(transform,amat) # sa(new) <- sa(old) <- ca

            ## multiple crystal poles may exist for each given (hkl)
            ## due to the crystal symmetry
            for j, pole_ca in enumerate(poles_ca):
                poles_projected[i,j,:] = np.dot(amat,pole_ca)
                poles_col[i,j,:]  = gr[4:] ## can be void for "TEX_PHx.OUT"
            poles_wgt[i,:]  = wgt

        poles_projected  = poles_projected.reshape( (len(grains)*len(poles_ca),3))
        poles_wgt = poles_wgt.reshape((len(grains)*len(poles_ca)))
        poles_col = poles_col.reshape((len(grains)*len(poles_ca),grains.shape[-1]-4))
    elif proj=='ipf':
        ## poles_projected can be either crystal poles (PF) or sample poles (IPF)

        # Now, in this case, pole is referring to sample direction.
        if csym=='cubic':
            H=sym.cubic()
        elif csym=='hexag':
            H=sym.hexag()
        elif csym=='ortho':
            H=sym.ortho()
        else:
            raise IOError('Not valid symmetry for pf')
        nsymop=H.shape[0]
        # empty np arrays
        poles_projected  = np.zeros((len(grains),nsymop*2,3)) #forward & backward
        poles_wgt        = np.zeros((len(grains),nsymop*2))
        poles_col        = np.zeros((len(grains),nsymop*2,grains.shape[-1]-4))
        for i, gr in enumerate(grains):
            phi1,phi,phi2,wgt=gr[:4]
            amat=euler(phi1,phi,phi2,a=None,echo=False) ## ca<-sa
            pole_ca=np.dot(amat, np.array(pole))## ca<-sa, sa
            for j, h in enumerate(H): # ca(new)<-ca(old)
                poles_projected[i,j*2,  :]=np.dot(h, pole_ca)
                poles_projected[i,j*2+1,:]=np.dot(h,-pole_ca)
                poles_col[i,j*2,:]=gr[4:]
                poles_col[i,j*2+1,:]=gr[4:]
            poles_wgt[i,:] = wgt

        #poles_projected[:,nsymop:2*nsymop,:]=-poles_projected[:,0:nsymop,:]
        #poles_col[:,nsymop:2*nsymop]=poles_col[:,0:nsymop]
        #poles_wgt[:,nsymop:2*nsymop]=poles_wgt[:,0:nsymop]
        ## reshaping
        poles_projected=poles_projected.reshape( (len(grains)*nsymop*2),3)
        poles_wgt=poles_wgt.reshape((len(grains)*nsymop*2))
        poles_col=poles_col.reshape((len(grains)*nsymop*2,grains.shape[-1]-4))


    if iopt==1:
        XY=[]
        WGT=[]
        COL_val=[]
        for ip, pole in enumerate(poles_projected):
            ## Convert each 3D pole to (x,y) coordinates
            x,y=projection(pole)
            if x**2+y**2<=1+tiny: ## If within the circle
                y=-y; x=-x
                XY.append([x,y])
                WGT.append(poles_wgt[ip])
                COL_val.append(poles_col[ip,:])

        return np.array(XY), np.array(WGT), np.array(COL_val)


    ## Full Sphere (-pi, +pi) and (0, pi)
    #x = np.arange(-180., 180.+tiny, dth) ## in-plane rotation
    #y = np.arange(   0., 180.+tiny, dph) ## tilting
    nx, ny = int(360./dth), int(180./dph)
    f = np.zeros((nx,ny))

    ## Semi Sphere (-pi, +pi) and (0, pi/2)
    x_node = np.arange(-180.,180.+tiny,dth) ## in-plane rotation
    y_node = np.arange(   0., 90.+tiny,dph) ## tilting
    nx_node = len(x_node); ny_node = len(y_node)
    nodes = np.zeros((nx_node,ny_node))
    f = pole2f(poles_projected,poles_wgt,dth,dph,f.copy())

    ## Normalization (m.u.r)
    fsum=f[:,:int(ny/2)].flatten().sum() ## sum of weights on the half-sphere
    z = np.zeros((ny+1))
    for i in range(ny+1): z[i] = np.pi/float(ny)*i
    dx_   = 2.*np.pi/nx
    dcosz = -np.diff(np.cos(z))
    fnorm = dcosz*dx_/(2*np.pi)
    f_ori=f.copy() # without normalization
    f     = f/fnorm/fsum

    ## Extension of f_bounds - see algorithm ipynb
    f_bounds = np.zeros((nx+2,ny+2))
    f_bounds[1:-1,1:-1]=f[ :, :]
    f_bounds[  -1,1:-1]=f[ 0, :]
    f_bounds[   0,1:-1]=f[-1, :]
    f_bounds[1:-1,   0]=f[ :,-1]
    f_bounds[   0,   0]=f_bounds[ 0,-2]
    f_bounds[  -1,   0]=f_bounds[-1,-2]
    f_bounds[   :,  -1]=f_bounds[ :, 1]

    ## pole-figure-weighted quantities for the extra columns in <esgr_x.out>
    ncols=grains.shape[-1]## addition to the pole weights
    if ncols>4:
        fcol=np.zeros((nx,ny,ncols-4))
        fcol = pole2f_cols(poles_projected,
                           poles_wgt,poles_col,f_ori,
                           dth,dph,fcol.copy())

        ## Extension of f_bounds - see algorithm ipynb
        f_bounds_col = np.zeros((nx+2,ny+2,ncols-4))
        for icol in range(ncols-4):
            f_bounds_col[1:-1,1:-1,icol]=fcol[ :, :,icol]
            f_bounds_col[  -1,1:-1,icol]=fcol[ 0, :,icol]
            f_bounds_col[   0,1:-1,icol]=fcol[-1, :,icol]
            f_bounds_col[1:-1,   0,icol]=fcol[ :,-1,icol]
            f_bounds_col[   0,   0,icol]=f_bounds_col[ 0,-2,icol]
            f_bounds_col[  -1,   0,icol]=f_bounds_col[-1,-2,icol]
            f_bounds_col[   :,  -1,icol]=f_bounds_col[ :, 1,icol]
        nodes_col = np.zeros((*nodes.shape,ncols-4))

    ## Use average of the four adjacent neighbouring nodes of pole figures
    for i in range(len(nodes)):
        for j in range(len(nodes[i])):
            nodes[i,j] = (f_bounds[i:i+2,j:j+2]).sum()/4.
            if ncols>4:
                for icol in range(ncols-4):
                    nodes_col[i,j,icol] = (f_bounds_col[i:i+2,j:j+2,icol]).sum()/4.

    ## Centeral region is using an avergage around the rim
    for i in range(n_rim):
        nodes[:,i]=(nodes[:,i].sum())/len(nodes[:,i])
        if ncols>4:
            for icol in range(ncols-4):
                nodes_col[:,i,icol]=(nodes_col[:,i,icol].sum())/len(nodes_col[:,i,icol])

    if ncols>4: # for those files in which columns in additon to (phi1, phi, phi2, wgt) are given
        return nodes, nodes_col
    else:
        return nodes

def __equiv__(miller=None, csym=None,
              cdim=[1.,1.,1.], cang=[90.,90.,90.]):
    """
    Provided the miller indices,
    Crystallographically equivalent and only unique
    vectors are returned.

    ---------
    Arguments
    ---------
    miller = None  , e.g. [1,1,1]
    csym   = 'cubic'
    cdim   = [ 1, 1, 1]
    cang   = [90,90,90]
    """
    start = time.time()
    from .sym import cv
    # from .sym import cubic, hexag, get_icsym, ortho
    #from sym_cy import cubic, hexag
    from . import sym    #python compiled
    #import sym_cy #cython compiled
    #from sym.py cvec, cubic, and hexgonal modules are brought in
    if type(miller)==type(None): raise IOError("Miller index should be given")

    vect = np.array(miller)
    norm = 0.; sneq = []
    temp = vect.copy()

    icsym=sym.get_icsym(csym)

    if csym=='cubic':
        H = sym.cubic()  #operators
        if False:
            for i in range(len(H)):
                sneq.append(np.dot(H[i], vect))
        else:
            v = cv(miller=vect,icsym=icsym,cdim=cdim,cang=cang)
            sneq=np.zeros((len(H),3))
            for m in range(len(H)):
                aux33=H[m].copy()
                bux3=np.zeros((3))
                for i in range(3):
                    for j in range(3):
                        bux3[i]=bux3[i]+aux33[i,j]*v[j]
                sneq[m,:]=bux3[:]
    elif csym=='hexag':
        H = sym.hexag() #operators
        v = cv(miller=vect, icsym=icsym,cdim=cdim, cang=cang)
        sneq=np.zeros((len(H),3))
        for m in range(len(H)):
            aux33=H[m].copy()
            bux3=np.zeros((3))
            for i in range(3):
                for j in range(3):
                    bux3[i]=bux3[i]+aux33[i,j]*v[j]
            sneq[m,:]=bux3[:]
    elif csym=='ortho':
        H = sym.ortho()  #operators
        v = cv(miller=vect, icsym=icsym,cdim=cdim, cang=cang)
        sneq=np.zeros((len(H),3))
        for m in range(len(H)):
            aux33=H[m].copy()
            bux3=np.zeros((3))
            for i in range(3):
                for j in range(3):
                    bux3[i]=bux3[i]+aux33[i,j]*v[j]
            sneq[m,:]=bux3[:]

    elif csym=='None':
        #H = [np.identity(3)]
        sneq = [vect]
    elif csym=='centro':
        sneq = [vect, -vect]
    else:
        print('Given symmetry, %s, is not expected'%csym)
        input('Enter to raise an error and quits the job');
        raise IOError

    #print 'elapsed time during v calculation: %8.6f'%
    #(time.time()-start)
    #####-------------------------------
    # start = time.time()
    stacked = [] #empty unique vectors
            # is cH in the already existing stacked list?
            # yes: pass
            # no : add

    ## filtering the sneq under whether or not it is unique
    for i in range(len(sneq)):
        cH = sneq[i].copy()  #current vector
        if __isunique__(a=cH, b=stacked):
            stacked.append(cH)

    ## if v[2] is minus, mutiply minus sign to the vector.
    for i in range(len(stacked)):
        if stacked[i][2]<0:
            stacked[i] = stacked[i]*-1
    #print 'elapsed time during the rest: %8.6f'%
    #(time.time()-start)
    return np.array(stacked)

### Excutes the module in the command line with arguments and options
def main(filename, pfmode, gr, csym):
    """
    plots the (100),(110) and (111) pole figures of cubic crystal.
    (to be including hexagonal)

    Arugments
    =========
      filename = '00500.cmb' for instance.
      csym : crystal symmetry 'cubic', 'hexag'
      pfmode ='contour','contourf', 'dot'
      """
    if gr!=None: a = polefigure(grains=gr, csym=csym)
    else: a = polefigure(filename=filename, csym=csym)
    #
    a.pf_new(mode=pfmode)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import getopt, sys

    ## arguments ------------------------------- ##
    try: opts, args = getopt.getopt(
        sys.argv[1:], 'm:i:o:c:s')#, 'output=', 'csym='])

    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    ## ----------------------------------------- ##
    ## default options
    ishow = False
    mode = 'contourf'
    csym = 'cubic'
    outputfile = 'temp.pdf'

    for o, a in opts:
        if o in ('-i'): inputfile = a
        elif o in ('-o'): outputfile = a
        elif o in ('-c'): csym = a
        elif o in ('-m'): mode = a
        elif o in ('-s'): ishow = True
        else: assert False, 'unhandled option'

    main(filename=inputfile, csym=csym, pfmode=mode)
    plt.gcf().savefig(outputfile)
    if ishow==True: plt.show()


def parse_epf(fn,n_unit=79):
    """
    Parse popLA-generated epf files

    Arguments
    ---------
    fn
    n_unit (the number of lines in an epf file)
    """
    nline_each_block = n_unit
    with open(fn,'r') as fo:
        string = fo.read()
        lines  = string.split('\n')[:-1]
        nl = len(lines)

    print('# of lines : %i'%nl)
    print('# of blocks: %i'%(float(nl)/float(n_unit)))

    nb = int(float(nl)/float(n_unit))

    blocks = []
    for i in range(nb):
        i0 = n_unit * i+1
        i1 = n_unit * (i+1)

        l = lines[i0:i1]

        l[0]=l[0][1:]
        # print '*'*10
        # print l[0]
        # print l[-1]

        b = ''
        # print 'number of lines:', len(l)
        for j in range(len(l)):
            b = '%s%s\n'%(b,l[j])
        blocks.append(b)

    return blocks

def axis2vect(i):
    """
    Argument
    -------
    <axis label> (1, 2, or 3; for the opposite directions use, -1, -2, or -3)

    Returns
    -------
    <Unit vector in the form of list array>
    """
    if i==1:
        return [1,0,0]
    if i==2:
        return [0,1,0]
    if i==3:
        return [0,0,1]
    if i==-1:
        return [-1,0,0]
    if i==-2:
        return [0,-1,0]
    if i==-3:
        return [0,0,-1]

def axes2transf(x,y):
    """
    Arguments
    ---------
    <x>:  (label of horizontal axis (right))
    <y>:  (label of vertical axis (top))

    ** note that, z (the direction made by x cross y) is naturally determined.

    Returns
    -------
    <mat> : transformation matrix that rotates given lab axes to
             the given (x, y, z)

          ** Use this transformation matrix to rotate the VPSC's discrete orientations (in Bunge)


    ## Example
    Refer to <8378_pf.ipynb> located under <matData/vpscData/AZ31_Wang_8378> of VPSCX repository.
    """
    xv=axis2vect(x)
    yv=axis2vect(y)
    zv=np.cross(xv,yv)
    mat=np.zeros((3,3)) ## xv<- e1, yv<- e2, zv<- e3
    mat[:,0]=xv[:]
    mat[:,1]=yv[:]
    mat[:,2]=zv[:]
    print('xv:',xv)
    print('yv:',yv)
    print('zv:',zv)
    return mat

def get_grid_from_angles(dth,dph,rot):
    """
    Calculate meshgrid used for pole figure contouring
    based on two angular increments (dth, dph). Additional
    in-plane rotation can be applied by using argument <rot>
    given in radian.

    Arguments
    ---------
    dth : in-plane rotation
    dph : tilting (half-sphere)
    rot : additional in-plane rotation (radian)

    Returns
    -------
    R   : meshgrid distance from the center (R) coordinates
    PHI : meshgrid rotation
    x   : meshgrid cartesian x (horizontal)
    y   : meshgrid cartesian y (vertical)
    """
    tiny = 1.e-9

    x_node = np.arange(-180.,180.+tiny,dth) ## in-plane rotation
    y_node = np.arange(   0., 90.+tiny,dph) ## tilting (half-sphere)
    XN, YN = np.meshgrid(x_node,y_node)

    #--------------------------------------------------#
    ## plotting / resolution
    nm     = int((360.0 - 0.)/dth) ## in-plane rotation
    nn     = int((180. - 90.)/dph) ## tilting
    theta  = np.linspace(pi, pi/2., nn+1)
    phi    = np.linspace(0., 2.*pi, nm+1)
    r      = np.sin(theta)/(1-np.cos(theta))
    R, PHI = np.meshgrid(r,phi)
    PHI    = PHI + rot ## default: rot=0.
    x      = R*np.cos(PHI); y = R*np.sin(PHI)

    return R, PHI, x, y


def xy_grid_from_shape(shape):
    """
    Obtain x,y grid of sphere on which (inverse) pole figures
    are drawn. Use 'shape' of gridded weight array, which might
    have been 'refined'.

    Argument
    --------
    shape

    Returns
    -------
    x,y
    """
    pi=np.pi
    nm,nn=np.array(shape,dtype='int')-1
    theta  = np.linspace(pi, pi/2., nn+1)
    phi    = np.linspace(0., 2.*pi, nm+1)
    r      = np.sin(theta)/(1-np.cos(theta))
    R, PHI = np.meshgrid(r,phi)
    x      = R*np.cos(PHI); y = R*np.sin(PHI)
    return x,y


def calc_arc(aca,rots):#a,b,fnsx,nang):
    """
    Arguments
    ---------
    aca crystal direction
    rots: rotation matrix

    Returns
    -------
    v_arc
    """
    nang=rots.shape[0]
    v_arc=np.zeros((nang,3))
    for i,rot in enumerate(rots):
        v_arc[i,:]=np.dot(rot,aca)
    return v_arc






def proj(a):
    """
    Apply stereographic projection of the given pole single pole 'a' in (3)
    or 'a' in (nvec,3) with 'nvec' being the number of separate poles

    Argument
    --------
    a : pole(s) (unit vector)

    Returns
    -------
    Multiple coordinates of 'nvec' number of projected poles (nvec,2)
        or coordinate of a single projected pole in (2,) shape
    """
    if len(a.shape)>1:
        vs=a.copy()
    elif len(a.shape)==1:
        vs=np.zeros((1,3))
        vs[0,:]=a[::]
    nvec=vs.shape[0]
    coords=np.zeros((nvec,2))
    for i, v in enumerate(vs):
        coords[i,:]=projection(v)
    if len(a.shape)>1: return coords
    elif len(a.shape)==1: return coords[0]


def get_ipf_boundary(nres=5,fnsx=None):
    """
    Given a, b, and c poles, calculate the fundamental triangle
    in the sphere in which the inverse pole figure contours are
    bounded. The algorithm is as follows:
    1. Pair up a,b, and c 'Miller-indexed' crystal plane normals such that
       (a,b), (b,c), (c,a)
    2. For each pair, calculates the


    Arguments
    ---------
    nres: the number of points belonging to each arc
          of (a,b), (b,c), and (c,a) pairs makes.
    fnsx: Name of single crystal file used in VPSC or dEVPSC code
        from which the crystallographic information is obtained.

    Returns
    -------
    coords: The coordinates of boundary used in the inverse pole figure.
    a, b, c : Miller-indexed poles which consist the fundamental triangle.
    """
    from . import sym

    #icsym=sym.get_icsym(csym)
    csym, cdim, cang = sym.read_fnsx(fnsx)

    # ## determine the three poles that define the fundamental zones.
    if csym=='cubic':
        a=[0,0,1]
        b=[1,0,1]
        c=[1,1,1]
    elif csym=='hexag':
        a=[0,0,0,2]
        b=[1,0,-1,0]
        c=[2,-1,-1,0]
    elif csym=='ortho':
        a=[0,0,1]
        b=[1,0,0]
        c=[0,1,0]
    else:
        raise IOError('Error: need to validate other crystal symmetries.')


    pairs=[[a,b],[b,c],[c,a]]
    coords=np.zeros((2,(nres-1)*3+1))

    for i, pair in enumerate(pairs[:3]):
        aca,bca,thf,vref,rots=sym.calc_vref_and_rot(*pair,fnsx,nres)
        varc=calc_arc(aca,rots)
        xy=proj(varc)
        i0=i*(nres-1)
        i1=i0+nres-1
        coords[:,i0:i1]=xy[0:-1,:].T
    coords[:,-1]=coords[:,0]
    return coords, a, b, c


def gen_mask_contour(boundary,x,y,shape):
    """
    Given the triangle boundary of inverse pole figure,
    determine coordinates that are NOT within the boundary.
    Then, return 'mask' in the shape of pole figure grid.
    The latter is to remove the domains outside the fundamental
    stereographic triangle.

    Arguments
    ---------
    boundary
    x
    y
    shape : shape of pole figure grid, which is used to contruct the 'mask' array

    Returns
    -------
    mask: masking array to removed the grid outside the boundary.
    """
    from shapely import geometry
    from shapely.geometry import Point, Polygon

    xb,yb=boundary
    poly=Polygon(zip(xb,yb))
    mask=np.empty(shape,dtype='object')
    for i in range(shape[0]):
        for j in range(shape[1]):
            p1,p2=x[i,j],y[i,j]
            if abs(p2)<1e-3: p2=1e-5 # a little trick to include the point on the horizontal line
            point=Point(p1,p2)
            mask[i,j]=not(point.within(poly))
    return mask

def get_within_triangle(boundary,XY):
    """
    Similar to <gen_mask_contour>, select XY points only within
    the given boundary.

    Arguments
    ---------
    boundary
    XY

    Returns
    -------
    XY_new
    tags
    """
    from shapely import geometry
    from shapely.geometry import Point, Polygon

    xb,yb=boundary
    poly=Polygon(zip(xb,yb))
    XY_new=[]
    tags=np.empty(len(XY),dtype='bool')
    for i, xy in enumerate(XY):
        x,y=xy
        point=Point(x,y)
        t=point.within(poly)
        tags[i]=t
        if t:
            XY_new.append(xy)
    return XY_new, tags



def get_pf_color_map_and_norm(levels,cmap,lev_norm_log, mns, mxs, nlev):
    """
    Based on the given set of arguments,
    obtain necessary objects for color-contouring.


    Arguments
    ---------
    levels,cmap,lev_norm_log, mns, mxs, nlev

    Returns
    -------
    levels, cm_norm, cmap_mpl, color_mapping
    """
    import numpy as np
    from matplotlib.colors import LogNorm
    import matplotlib.cm
    import matplotlib.pyplot

    if type(levels)==type(None):
        if lev_norm_log:
            ## To prevent log (0) -> np.nan
            ## hardwire minimum value
            if mns==0: mns = 0.5
            levels = np.logspace(
                np.log10(mns),np.log10(mxs),nlev)
            cm_norm = LogNorm()
        else:
            levels = np.linspace(mns,mxs,nlev)
            cm_norm = None
    else:
        cm_norm = None

    try:
        cmap_mpl = matplotlib.cm.get_cmap(cmap)
    except:
        # print("**Warning: Couldn't use matplotlib.cm.get_cmap")
        # print(' I am now using matplotlib.pyplot.get_cmap')
        cmap_mpl = matplotlib.pyplot.get_cmap(cmap)

    color_mapping = matplotlib.cm.ScalarMappable(
        norm=cm_norm,cmap=cmap_mpl)

    return  levels, cm_norm, cmap_mpl, color_mapping
