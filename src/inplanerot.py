def inplane_rotate_fntex(fntx,ang):
    """
    Arguments
    ---------
    fntx
    ang

    Returns
    -------
    fntex_rotated: file in which rotated texture is stored
    grains_rot:    numpy array of the discrete orientions
    """
    from .euler import euler
    import numpy as np
    import tempfile

    R=np.zeros((3,3)) ## new_sa <- old_sa
    cos=np.cos
    R[0,0]=cos(np.deg2rad(ang))
    R[0,1]=cos(np.deg2rad(90.-ang))
    R[1,0]=cos(np.deg2rad(90.+ang))
    R[1,1]=cos(np.deg2rad(ang))
    R[2,2]=1.
    R=R.T ## old_sa <- new_sa

    grains=np.loadtxt(fntx,skiprows=4).T
    grains_rot=np.zeros(grains.shape)

    for i, gr in enumerate(grains.T):
        phi1,phi,phi2,wgt=gr[:]
        g=euler(ph=phi1,th=phi,tm=phi2,echo=False) # ca <- sa
        newg=np.dot(g,R)
        phi1,phi,phi2=euler(a=newg,echo=False)
        grains_rot[:,i]=phi1,phi,phi2,wgt

    fntex_rotated=tempfile.mktemp()

    with open(fntex_rotated,'w') as fo:
        fo.write(f'dum\ndum\ndum\nB {grains.shape[-1]}\n')
        for i, gr in enumerate(grains_rot.T):
            fo.write('%9.3f %9.3f %9.3f %9.5f\n'%(gr[0],gr[1],gr[2],gr[3]))

    return fntex_rotated, grains_rot
