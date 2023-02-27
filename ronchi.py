import numpy as np
import scipy.interpolate as si
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('qtagg')
plt.style.use('dark_background')

def ronchi(
    n: int,
    m: int,
    d: float,
    phi: float,
    xoff: float,
    neg: bool=False
):
    """Generate a ronchi ruling.

    Parameters
    ----------
    n : int
        Size along x of the ronchi ruling, in pixels
    m : int
        Size along y of the ronchi ruling, in pixels
    d : float
        Width of the ronchi ruling bands, in pixels
    phi : float
        Orientation angle of the ronchi grating
    xoff : float
        ???
    neg : bool
        If true, invert the ronchi ruling

    Returns
    -------


    """
    x, y = np.meshgrid(np.arange(n), np.arange(m))
    xp = x*np.cos(phi*np.pi/180) + y*np.sin(phi*np.pi/180) + xoff
    b = np.floor(xp/d)
    return np.bitwise_xor(~np.abs(b % 2).astype(bool), neg)

def ronchigram(a, f, d, zoff, phi, sph, coma, astig, defocus, n, m, neg=False, xsize=1, xoff=1):
    """Generate a ronchigram.

    Parameters
    ----------
    a : float
        Radius of the exit pupil
    f : float
        Focal length of system under test
    d : int
        Width of the ronchi ruling bands, in pixels
    zoff : float
        Offset from the focal point along the optical axis
    phi : float
        Orientation angle of the ronchi grating
    sph : float
        Spherical aberration coefficient
    coma : float
        Coma aberration coefficient
    astig : float
        Astigmatism aberration coefficient
    defocus : float
        Defocus coefficient
    n : int
        Size along x of the ronchi ruling, in pixels
    m : int
        Size along y of the ronchi ruling, in pixels
    neg : bool
        If true, invert the ronchi ruling
    xsize : float
        ???
    xoff : float
        ???
    """
    # import pudb; pudb.set_trace()
    # Normalize variables
    fn = f/a
    zn = zoff/a
    dn = d/a

    # Calculate defocus

    # Calculate radius of optical system
    rn = 2*fn

    # Calculate magnification; this is needed because the ronchi ruling is to be
    # defined at the exit pupil of the system.
    if (zn == 0):
        mag = -rn
    else:
        mag = -rn/zn

    # define the normal coordinates
    x, y = np.meshgrid(
        np.linspace(-1, 1, n),
        np.linspace(-1, 1, m)
    )

    # calculate the aberrated coordinates (i.e., the difference in path length)
    xa = x-mag*rn*(4*sph*x*(x**2 + y**2) + 2*coma*x*y + 2*astig*x + 2*defocus*x)
    ya = y-mag*rn*(4*sph*x*(x**2 + y**2) + coma*x**2 + 3*coma*y**2 + 6*astig*y + 2*defocus*y)

    # define the ronchi ruling pitch in pixels
    d1max = np.max(xsize*dn*np.abs(mag))
    d1 = np.round(np.where(d1max > 1, d1max, 1))

    # set a maximum normalized aberrated coordinates big enough to accommodate the abberated coordinates
    c = 2

    # redefine the magnified coordinates
    xm, ym = np.meshgrid(
        np.linspace(-c, c, c*n),
        np.linspace(-c, c, c*m)
    )

    # xoff in pixels
    o1 = c*xsize*mag*xoff

    # calculate the magnified ronchi ruling
    g1 = ronchi(
        n=c*n,
        m=c*m,
        d=d1,
        phi=phi,
        xoff=o1,
        neg=neg
    )

    # interpolate the ruling function on aberrated coordinates
    return si.griddata(
        (xm.flatten(), ym.flatten()),
        g1.flatten(),
        (xa, ya)
    )


fig, ax = plt.subplots(1, 1, figsize=(10, 10))

nxy = 200
rg = ronchigram(
    a=25,
    f=200,
    d=0.1,
    zoff=-30,
    phi=0,
    sph=0,
    coma=0,
    astig=0,
    defocus=0.0005,
    n=nxy,
    m=nxy,
    neg=False,
    xsize=1,
    xoff=0,
)

ax.imshow(rg, interpolation='none', origin='lower')
plt.show()
