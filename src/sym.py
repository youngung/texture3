"""
General crystallographic symmetry operators
Cubic and hexagonal operators are available.


-- List of symmery operators
def cubic()
def hexag()
def tetra()
def triclinic()
"""

import numpy as np
import time
## symmetry operators





def __60_120_rot111__(h):
    """
    For the given h operation,
    rotations of (pi/3) & (2*pi/3) around <111>
    are performed and returned

    *cubic
    """
    hx = h.copy()
    h60 = np.zeros((3,3)); h120 = np.zeros((3,3))
    h60[0,2] = 1.
    h60[1,0] = 1.
    h60[2,1] = 1.

    h120[0,1] = 1.
    h120[1,2] = 1.
    h120[2,0] = 1.
    return np.dot(h60,hx), np.dot(h120,hx)

def __mirror_110__(h):
    """
    Given the operation h, mirrored across the (110) plane returned

    *cubic
    """
    hx = h.copy()
    hm = np.zeros((3,3))
    hm[0,1] = 1.
    hm[1,0] = 1.
    hm[2,2] = 1.
    return np.dot(hm, hx)

def __rot_90_180_270__(h):
    """
    Given the operation h,
    the three rotated operations are returned

    *cubic
    """
    cos = np.cos; sin = np.sin; pi = np.pi
    hx = np.zeros((3,3,3))
    h_ = h.copy(); htemp = []

    for m in range(3):
        ang = pi/2. * float(m+1)
        hx[m,0,0] = cos(ang)
        hx[m,1,1] = cos(ang)
        hx[m,2,2] = 1.0
        hx[m,0,1] = -sin(ang)
        hx[m,1,0] = sin(ang)
        hx[m,0,2] = 0.
        hx[m,2,0] = 0.
        hx[m,1,2] = 0.
        hx[m,2,1] = 0.
        pass
    for m in range(3):
        htemp.append( np.dot(hx[m], h_) )
        pass
    return np.array(htemp)

def __rot_nrot_x1__(h,nrot):
    """
    Mirror plane at 30 or 60 or 45 deg with respect to x1

    *hexagonal, trigonal, tetragonal

    hexa: nrot = 6
    trig: nrot = 3
    tetr: nrot = 4
    """
    cos = np.cos; sin = np.sin; pi=np.pi
    hx = np.zeros((3,3))
    ang = pi/float(nrot)
    hx[0,0] = cos(ang)**2 - sin(ang)**2
    hx[1,1] = -hx[0,0]
    hx[2,2] = 1.0
    hx[0,1] = 2.*cos(ang)*sin(ang)
    hx[1,0] = hx[0,1]

    #print('rot_nrot_x1')
    #for j in range(3):
    #print('%5.2f %5.2f %5.2f'%(hx[j][0],hx[j][1],hx[j][2]))
    #print('--')

    return np.dot(hx,h)

def __rot_nrot_001__(h, csym=None):
    """
    Rotations of 2*pi/nrot around axis <001>

    *hexagoanl, trigonal, tetragonal

    ---------
    Arguments
    h: symmetry operators
    csym: 'hexa'
    """
    if   csym=='hexag': nrot=6
    elif csym=='trigo': nrot=3
    elif csym=='tetra': nrot=4
    else: print('Unexpected Error'); raise IOError

    cos = np.cos; sin = np.sin; pi = np.pi
    hx = np.zeros((nrot-1,3,3))
    h_ = h.copy(); htemp = []

    for nr in range(nrot-1):
        ang = (nr+1)*2.*pi/nrot
        hx[nr,0,0] = cos(ang)
        hx[nr,1,1] = cos(ang)
        hx[nr,2,2] = 1.0
        hx[nr,0,1] =-sin(ang)
        hx[nr,1,0] = sin(ang)

    print('nx for nrot-1',nrot-1)
    for nr in range(nrot-1):
        for i in range(3):
            print('%5.2f %5.2f %5.2f'%(hx[nr,i,0],hx[nr,i,1],hx[nr,i,2]))
        print('--')
    print()
    print()

    for nr in range(nrot-1):
        htemp.append(np.dot(hx[nr], h_))

    print('htemp:')
    for nr in range(nrot-1):
        for i in range(3):
            print('%5.2f %5.2f %5.2f'%(htemp[nr][i,0],htemp[nr][i,1],htemp[nr][i,2]))
        print('--')



    return np.array(htemp)

def __trim0__(h):
    """
    if a value in the matrix h is fairly close to +-0.
    then returns zero. In that way, trimming is performed
    on the every component of given h matrix.
    """
    hx = h.copy()
    for i in range(len(hx)):
        for j in range(len(hx[i])):
            if abs(hx[i,j]) < 1e-6:
                hx[i,j] = 0.
    return hx


""" -- mmm sample symmetry is found in COD_conv.py """
def __mmm__():
    m0 = [[ 1, 0, 0], [0, 1, 0], [0, 0, 1]]
    m1 = [[ 1, 0, 0], [0,-1, 0], [0, 0,-1]]
    m2 = [[-1, 0, 0], [0, 1, 0], [0, 0,-1]]
    m3 = [[-1, 0, 0], [0,-1, 0], [0, 0, 1]]
    h = np.array([m0,m1,m2,m3])
    return h

### deprecated ###
# def __ortho__(v):
#     """
#     Orthogonal sample symmetry to a vector (v) in 3D
#     """
#     v1 = np.aray([a[0], a[1], a[2]])
#     v2 = np.array([-a[0], a[1], a[2]])
#     v3 = np.array([a[0], -a[1], a[2]])
#     v4 = np.array([-a[0], -a[1], a[2]])

#     v5 = v1.copy()* -1
#     v6 = v2.copy()* -1
#     v7 = v3.copy()* -1
#     v8 = v4.copy()* -1

#     return v1, v2, v3, v4


## cubic symmetry
def cubic():
    H = []     # H is the master list containing all the numpy arrays of operations
    H.append(np.identity(3))    # identity operation

    # rotations of (pi/3) & (2*pi/3) around <111>
    niter = len(H)
    for i in range(niter):
        h60, h120 = __60_120_rot111__(h=H[i].copy())
        h0 = h60.copy(); h1 = h120.copy()
        H.append(h0)
        H.append(h1)

    # mirror across the plane (110)
    niter = len(H)
    for i in range(niter):
        h = __mirror_110__(h=H[i].copy())
        H.append(h)

    # rotations of 90, 180, 270 around x3
    niter = len(H)
    for i in range(niter):
        h1, h2, h3 = __rot_90_180_270__(h=H[i].copy())
        h90 = h1.copy(); h180 = h2.copy(); h270 = h3.copy()
        H.append(h90)
        H.append(h180)
        H.append(h270)

    H = np.array(H) # Make the H as numpy array

    # trim the values.
    for i in range(len(H)):
        H[i] = __trim0__(h=H[i])
    return H

def cubic_centro():
    h_old = cubic()
    h_new = []
    h_n = [[-1,0,0],[0,-1,0],[0,0,-1]]
    for i in range(len(h_old)):
        h_new.append(np.dot(h_old[i],h_n))
    return h_new


def triclinic():
    H = []
    H.append(np.identity(3))
    return H

## hexagonal
def hexag():
    H = []
    H.append(np.identity(3))

    #mirror plane at 30 degree with respect to x1
    nrot = 6
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_x1__(h=H[i].copy(),nrot=nrot)
        H.append(h)

    hx=np.zeros((3,3,24))

    for i in range(nrot-1):
        nr=i+1
        ang=nr*2*np.pi/nrot
        # print('ang:',ang)
        hx[0,0,i]=np.cos(ang)
        hx[1,1,i]=np.cos(ang)
        hx[2,2,i]=1.
        hx[0,1,i]=-np.sin(ang)
        hx[1,0,i]= np.sin(ang)

    HS=[]
    for i in range(nrot-1):
        a=hx[:,:,i].copy()
        for k in range(len(H)):
            b=H[k].copy()
            aux=np.zeros((3,3))
            for m in range(3):
                for n in range(3):
                    for o in range(3):
                        aux[m,n]=aux[m,n]+a[m,o]*b[o,n]

            HS.append(aux)

    for i in range(len(HS)):
        H.append(HS[i])

    return np.array(H)

## orthorhombic
def ortho():
    H=[]
    H.append(np.identity(3))
    niter=len(H)

    pi=np.pi
    cp=np.cos(pi)
    sp=np.sin(pi)
    # 180 deg rotation around (001)
    h=np.zeros((3,3))
    h[0,0]=cp
    h[1,1]=cp
    h[2,2]=1.
    h[0,1]=-sp
    h[1,0]=sp
    H.append(h)

    # x-mirror & y-mirror
    h=np.identity(3)
    h[0,0]=-1
    H.append(h)
    h=np.identity(3)
    h[1,1]=-1
    H.append(h)
    return np.array(H)

## trigonal
def trigo():
    H = []
    H.append(np.identity(3))
    #mirror plane 60 degree with respect to x1
    nrot = 3
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_x1__(h=H[i].copy(), nrot=3)
        H.append(h)
    #rotations of 2*pi/3 around axis <001> for trigonals
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_001__(h=H[i], csym='trigo')
        H.append(h)

    for i in range(len(H)):
        H[i] = __trim0__(h=H[i])
    return H

## tetragonal
def tetra():
    H = []
    H.append(np.identity(3))
    #mirror plane at 45 degree with respect to x1
    nrot = 4
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_x1__(h=H[i].copy(), nrot=nrot)
        H.append(h)

    #rotations of 2*pi/4 around axis <001> for hexagonals.
    niter = len(H)
    for i in range(niter):
        h = __rot_nrot_001__(h=H[i], csym='tetra')
        for ix in range(len(h)):
            H.append(h[ix])

    for i in range(len(H)):
        H[i] = __trim0__(h=H[i])
    return H

##


def cvec(cdim=None, cang=None):
    """
    Generates and returns 'cvec[i,n]' of the unit cell that
    is characterized by unit cell dimension together with
    axes' angles

    ---------
    Arguments
      cdim=[1.,1.,1.]
      cang=[90.,90.,90.] : Should be in angle [90.,90.,90.] not radian
    """

    cdim = np.array(cdim)
    cang = np.array(cang)
    # angle to radian
    cang = cang * np.pi/180.
    # cvec
    cvec = np.zeros((3,3))

    cvec[0,0] = np.sin(cang[1])
    cvec[1,0] = 0.
    cvec[2,0] = np.cos(cang[1])

    cvec[0,1] = (np.cos(cang[2])-np.cos(cang[0])\
                     *np.cos(cang[1]))/np.sin(cang[1])
    cvec[2,1] = np.cos(cang[0])
    cvec[1,1] = np.sqrt(1.-cvec[0,1]**2-cvec[2,1]**2)

    cvec[0,2] = 0.
    cvec[1,2] = 0.
    cvec[2,2] = 1.

    # print('cvec')
    # for i in range(3):
    #     print('%5.2f %5.2f %5.2f'%(cvec[i,0],cvec[i,1],cvec[i,2]))
    # print('--')

    for i in range(3):
        for j in range(3):
            cvec[i,j] = cdim[j] * cvec[i,j]


    # print('cvec = cdim * cvec')
    # for i in range(3):
    #     print('%5.2f %5.2f %5.2f'%(cvec[i,0],cvec[i,1],cvec[i,2]))
    # print('--')

    return cvec

def get_icsym(crysym):
    """
    Covert cysym to icrysym to be compatible with 'crystal_symmetry.f' syntax.

    Argument
    --------
    crysym

    Returns
    -------
    icrysym
    """
    icrysym=0
    if crysym[:5].lower()=='cubic': icrysym=1
    if crysym[:5].lower()=='hexag': icrysym=2
    if crysym[:5].lower()=='trigo': icrysym=3
    if crysym[:5].lower()=='tetra': icrysym=4
    if crysym[:5].lower()=='ortho': icrysym=5
    if crysym[:5].lower()=='monoc': icrysym=6
    if crysym[:5].lower()=='tricl': icrysym=7


    if icrysym==0:
        raise IOError(f'Unexpected crysym is given {crysym}')
    return icrysym

def cv(miller, icsym=None, cdim=None, cang=None):
    """
    Creates a vector of the (plane normal) pole taking care of its unit cell's
    dimension and axes' angles.

    Arguments
    ---------
    pole : miller indexed pole
    icsym: 1 (cubic), 2 (hexag), 3 (trigo), 4 (tetra), 5 (ortho), 6 (monoc), 7 (tricl)

    Returns
    -------
    Cartesian vector corresponding to the given plane normal (miller).
    It is assumed that the plane normal is given in terms of Miller-indexed form.
    """
    if type(icsym)==type(None):
        raise IOError('icsym should be given to <sym.cv>')

    pole=miller.copy()

    if icsym==2 or icsym==3:
        pole[2]=pole[3]

        # below is for 'directions' not for plane normals.
        # pole[0]=pole[0]-pole[2]
        # pole[1]=pole[1]-pole[2]
        # pole[2]=pole[3]

    sqrt = np.sqrt
    cvect = cvec(cdim=cdim, cang=cang)

    s = np.zeros((3,))
    s[2] = ( pole[2]                                     ) / cvect[2,2]
    s[0] = ( pole[0] - cvect[2,0]*s[2]                   ) / cvect[0,0]
    s[1] = ( pole[1] - cvect[0,1]*s[0] - cvect[2,1]*s[2] ) / cvect[1,1]

    norm = sqrt(s[0]**2 + s[1]**2 + s[2]**2)
    for i in range(3):
        s[i] = s[i]/norm
        if abs(s[i])<1e-6: s[i]=0.
    return s


def read_fnsx(fnsx):
    """
    Read the single crystal file, and return crystal symmetry
    including other dimensions such as 'cang' and 'cdim'.

    Argument
    --------
    fnsx

    Returns
    -------
    csym
    cang
    cdim
    """
    with open(fnsx,'r',errors="ignore") as fo:
        lines=fo.read().split('\n')[:20]

    csym=lines[1].split()[0][:5].lower()
    cdim=np.array(lines[2].split()[:3],dtype='float')
    cang=np.array(lines[2].split()[3:6],dtype='float')
    return csym, cdim, cang

def calc_cvec(miller,fnsx):
    """
    Get a crystal vector of a Miller-indexed "plane" normal

    Arguments
    ---------
    miller: Miller-indexed crystal plane
    fnsx: Name of single crystal file used in VPSC or dEVPSC code
        from which the crystallographic information is obtained.

    Returns
    -------
    The normal vector of the given crystal plane (hkl).
    """
    # from TX import sym
    csym,cdim,cang = read_fnsx(fnsx)
    icsym=get_icsym(csym)
    _p_=cv(miller,icsym,cdim,cang)
    return _p_



def calc_vref_and_rot(a,b,fnsx,nang):
    """
    Given a miller-indexed crystal plane normal (a and b),
    calclulate the v ref and rotation matrices vref.

    This is used mainly to obtain the great circle that passes through
    the given a and b poles that are in Miller index. The crystal symmetry
    information from <fnsx> is used to construct the cartesian vectors
    out of the two miller-indexed poles (i.e., a and b).

    Arguments
    ---------
    a      : miller-indexed plane normal (hkl)
    b
    fnsx
    nang

    Returns
    -------
    aca   : cartesian vector of the given Miller-indexed 'a' pole (plane-normal)
    bca   : cartesian vector of the given Miller-indexed 'b' pole (plane-normal)
    thf   : the maximum rotation angle that connects <aca> to <bca>
    vref  : the axis about which the rotation occurs to connect aca and bca via the great circle.
    rots  : rotation matrices in the shape of (nang, 3, 3). 'nang' is the number of
            increments of the total rotation that connects aca to bca along the great circle.
    """
    from . import bcc_rolling_fiber
    aca=-calc_cvec(miller=a,fnsx=fnsx)
    bca=-calc_cvec(miller=b,fnsx=fnsx)
    vref=np.cross(aca,bca)

    ##
    thf=np.arccos(np.dot(aca,bca))
    ths=np.linspace(0,thf,nang)
    rots=np.zeros((nang,3,3))
    varc=np.zeros((nang,3))

    for i, th in enumerate(ths):
        rots[i,:,:]=bcc_rolling_fiber.vector_ang(vref,np.rad2deg(th))
    return aca, bca, thf, vref, rots
