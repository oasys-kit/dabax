

import numpy
from dabax.dabax_xraylib import DabaxXraylib
import xraylib

if __name__ == "__main__":
    dx = DabaxXraylib(file_f1f2='f1f2_Henke.dat')


    #
    # refractive index
    #
    dens = dx.element_density("Be")
    print("Refractive_Index_Re Be  dabax: ",  dx.Refractive_Index_Re ("Be",18.0, dens))
    print("Refractive_Index_Im Be  dabax: ",  dx.Refractive_Index_Im ("Be",18.0, dens))

    # loops
    energies = numpy.linspace(5,18,100)
    r1 = dx.Refractive_Index_Re("Be", energies, dens)
    i1 = dx.Refractive_Index_Im("Be", energies, dens)
    xr1 = numpy.zeros_like(energies)
    xi1 = numpy.zeros_like(energies)
    for i,energy in enumerate(energies):
        print("   delta @ %g keV, dabax: %g" % (energy, 1-r1[i]))
        xr1[i] = xraylib.Refractive_Index_Re ("Be", energies[i], dens)
        xi1[i] = xraylib.Refractive_Index_Im ("Be", energies[i], dens)

    from srxraylib.plot.gol import plot
    plot(energies, r1,
         energies, xr1,
         title="Real part of refraction index for Be",
         legend=["DABAX (Henke data)", "XRAYLIB"])
    plot(energies, i1,
         energies, xi1,
         title="Imaginary part of refraction index for Be",
         legend=["DABAX", "XRAYLIB"])
