# See also the examples are in the main of
#
#  dabax_base.py
#

import numpy
from dabax.dabax_base import DabaxBase

if __name__ == "__main__":
    dx = DabaxBase()

    #
    # misc
    #
    if True:
        from dabax.common_tools import atomic_symbols, atomic_names, atomic_number

        print(atomic_symbols()[14], atomic_names()[14])
        #
        # misc
        #
        print("Ge, Si: ", dx.atomic_weights(["Ge","Si"]))
        print("70Ge: ", dx.atomic_weights("70Ge"))

        print("Si atomic mass", dx.atomic_constants("Si", return_item=2))
        print("Si,Ge atomic mass", dx.atomic_constants(["Si", "Ge"], return_item=2))
        print("Si,Co atomic mass", dx.atomic_constants(["Si", "Co"], return_label='AtomicMass'))

        print("Z=27", atomic_symbols()[27])
        print("Ge Z=%d" % atomic_number("Ge"))

        print("Density Si: ", dx.element_density("Si"))

    #
    # f1f2
    #
    if True:
        energy, f1, f2 = dx.f1f2_extract("Si")

        energy_i = numpy.linspace(10,15000,200)
        f1_i, f2_i = dx.f1f2_interpolate("Si", energy=energy_i)
        print(">>>>", energy.shape, f1.shape, f2.shape)
        from srxraylib.plot.gol import plot
        plot(energy, f1,
             energy, f2,
             energy_i, f1_i,
             energy_i, f2_i,
             xlog=True, ylog=True, title="f1f2 Si",
             legend=['f1','f2','f1_i','f2_i'],
             marker=[None,None,'+','+'],
             linestyle=[None,None,'',''])




    #
    # CrossSec
    #
    if True:
        energy, cs = dx.crosssec_extract("Si")

        energy_i = numpy.linspace(10,15000,200)
        cs_i = dx.crosssec_interpolate("Si", energy=energy_i)
        print(">>>>", energy.shape, cs.shape)
        from srxraylib.plot.gol import plot
        plot(energy, cs,
             energy_i, cs_i,
             xlog=True, ylog=True, title="crosssec Si",
             legend=['cs','cs_i'],
             marker=[None,'+'],
             linestyle=[None,''])

        print(">>>> cs of Si at 10 kev %g barn/atom ", dx.crosssec_interpolate("Si", energy=50000.0))