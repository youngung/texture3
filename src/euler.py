import numpy as np
import math

def in_plane_rot(th):
    """
    Return an in-plane rotation matrix
    | cos(th)   -sin(th)   0 |
    | sin(th)    cos(th)   0 |
    | 0          0         1 |

    Argument
    --------
    th  (rotation angle in degree)
    """
    th = th * np.pi / 180.

    sth = np.sin(th)
    cth = np.cos(th)
    rot = np.zeros((3,3))
    rot[0,0] =  cth
    rot[0,1] = -sth
    rot[1,0] =  sth
    rot[1,1] =  cth
    rot[2,2] =  1.
    return rot

def euler(ph=None, th=None, tm=None, a=None, echo=True):
    """
    note:
          This is a pythonized fortran subroutine embedded in the VPSC7.sub
          done by youngung 2010-12-23

          Matrix A(i,j) is returned as a numpy array
          all Euler angles are in degree not in radian

          # A[i,j] transforms from SA to CA

          Thanks to python's non-hardwirable-arugment feature,
          if a matrix is not given, it automatically calculates
          the matrix with given angles.
          Vice versa, if matrix is given, given angle aurgments
          are ignored and new euler angles are returned.

    Nomenclature of Euler angle follows Bunge's convention
          ph = phi1,
          th = phi
          tm = phi2
    """
    if type(a).__name__=='NoneType':  a=np.resize(np.array(()),(3,3));iopt=2
    else:
        if type(a).__name__=='ndarray':
            iopt = 1
            pass
        else:
            print('Error: Unexpected Matrix a type')
            print('It should be numpy.ndarry!')
            raise IOError

    if iopt==1 :
        th = np.arccos(a[2,2])  #Radian
        if abs(a[2,2]) > 0.99999:
            tm = 0.
            ph = np.arctan2(a[0,1],a[0,0]) #Radian
        else:
            sth = np.sin(th)
            tm = np.arctan2(a[0,2]/sth,a[1,2]/sth)
            ph = np.arctan2(a[2,0]/sth,-a[2,1]/sth)
        th = th * 180./np.pi
        ph = ph * 180./np.pi
        tm = tm * 180./np.pi
        return [ph,th,tm] #phi1, phi, phi2

    elif (iopt == 2):
        angles = [ph,th,tm]
        if any(angles[i] == None for i in range(len(angles))):
            print('Angles must be give if iopt==2')
            raise IOError

        """ Convert the angle into Radian"""
        ph = ph * np.pi / 180.
        th = th * np.pi / 180.
        tm = tm * np.pi / 180.

        sph = np.sin(ph)
        cph = np.cos(ph)
        sth = np.sin(th)
        cth = np.cos(th)
        stm = np.sin(tm)
        ctm = np.cos(tm)

        a[0,0] =  ctm * cph - sph * stm * cth
        a[1,0] = -stm * cph - sph * ctm * cth
        a[2,0] =  sph * sth
        a[0,1] =  ctm * sph + cph * stm * cth
        a[1,1] = -sph * stm + cph * ctm * cth
        a[2,1] = -sth * cph
        a[0,2] =  sth * stm
        a[1,2] =  ctm * sth
        a[2,2] =  cth

        if echo==True:
            print('Matrix a is ')
            print(a)
        return a

    else: print('Error: Improper iopt input'); raise IOError



def eulers(phs=None, ths=None, tms=None, amats=None, echo=True,iopt=None):
    """
    This module is to allow calculate Euler angles or transformation matrices
    multiple times without requiring going over loops. The procedure relies
    on NumPy ufunc
    """
    # if type(a).__name__=='NoneType':  a=np.resize(np.array(()),(3,3));iopt=2
    # else:
    #     if type(a).__name__=='ndarray':
    #         iopt = 1
    #         pass
    #     else:
    #         print('Error: Unexpected Matrix a type')
    #         print('It should be numpy.ndarry!')
    #         raise IOError

    if type(iopt)==type(None):
        print('** Error: iopt should be given to eulers.')
        raise IOError('** Error')

    if iopt==1 :
        ths = np.arccos(amats[2,2])  #Radian
        tiny=1e-6
        flg=abs(amats[:,2,2]) > 1-tiny
        tms[flg,:]=0.
        phs[flg,:]=np.arctan2(amats[:,0,1],amats[:,0,0])
        # if abs(a[2,2] > 0.99999):
        #     tm = 0.
        #     ph = np.arctan2(a[0,1],a[0,0]) #Radian
        #else:
        sth = np.sin(ths)
        tms = np.arctan2(a[:,0,2]/sth,a[:,1,2]/sth)
        phs = np.arctan2(a[:,2,0]/sth,-a[:,2,1]/sth)
        ths = ths * 180./np.pi
        phs = phs * 180./np.pi
        tms = tms * 180./np.pi
        return [phs,ths,tms] #phi1, phi, phi2

    elif (iopt == 2):
        # angles = [phs,ths,tms]
        # if any(angles[i] == None for i in range(len(angles))):
        #     print('Angles must be give if iopt==2')
        #     raise IOError

        """ Convert the angle into Radian"""
        amats=np.zeros((len(phs),3,3))

        print(f'amats.shape: {amats.shape}')

        phs = phs * np.pi / 180.
        ths = ths * np.pi / 180.
        tms = tms * np.pi / 180.

        sph = np.sin(phs)
        cph = np.cos(phs)
        sth = np.sin(ths)
        cth = np.cos(ths)
        stm = np.sin(tms)
        ctm = np.cos(tms)

        amats[:,0,0] =  ctm * cph - sph * stm * cth
        amats[:,1,0] = -stm * cph - sph * ctm * cth
        amats[:,2,0] =  sph * sth
        amats[:,0,1] =  ctm * sph + cph * stm * cth
        amats[:,1,1] = -sph * stm + cph * ctm * cth
        amats[:,2,1] = -sth * cph
        amats[:,0,2] =  sth * stm
        amats[:,1,2] =  ctm * sth
        amats[:,2,2] =  cth

        if echo==True:
            print('Matrix a[0,:,:] is ')
            print(amats[0,:,:])
        return amats

    else: print('Error: Improper iopt input'); raise IOError
