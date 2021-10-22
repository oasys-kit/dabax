import numpy
from dabax.dabax_util import get_dabax_file
from dabax.dabax_util import get_f0_coeffs_from_dabax_file, calculate_f0_from_f0coeff, f0_with_fractional_charge
from dabax.dabax_util import f1f2_calc_dabax, atomic_weights_dabax, atomic_symbols_dabax, atomic_names_dabax
from dabax.dabax_util import atomic_constants_dabax, atomic_number_dabax, element_density_dabax
from dabax.dabax_util import cross_calc_dabax

from dabax.dabax_util import Crystal_GetCrystalsList, Crystal_GetCrystal, Crystal_dSpacing, Bragg_angle, CompoundParser

if __name__ == "__main__":

    import xraylib # for comparisons
    from srxraylib.plot.gol import plot

    # redefine the default server at ESRF because default server name is different outside and inside ESRF
    import socket
    if socket.getfqdn().find("esrf") >= 0:
        # dabax_repository = "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"
        dabax_repository = "/scisoft/DABAX/data"


    #
    # crystal tests
    #
    if True:
        print(get_dabax_file("Crystals.dat", dabax_repository=dabax_repository, verbose=0))

        print(get_f0_coeffs_from_dabax_file("Y3+",
                                            filename="f0_InterTables.dat",
                                            dabax_repository=dabax_repository))

        print(Crystal_GetCrystalsList(dabax_repository=dabax_repository))

        yb = Crystal_GetCrystal('YB66', filename='Crystals.dat', dabax_repository=dabax_repository)

        si = Crystal_GetCrystal("Si", dabax_repository=dabax_repository)
        print("Si 111 d-spacing: ", Crystal_dSpacing(si,1,1,1))
        print("Si 111 bragg angle at 10 keV [deg]: ", 180 / numpy.pi * Bragg_angle(si,10, 1,1,1))

    #
    # crystal vs xraylib tests
    #

    if True:
        print(Crystal_GetCrystal(entry_name='YB66', filename='Crystals.dat', dabax_repository=dabax_repository))

        # compare with xraylib
        xdabax = Crystal_GetCrystal(entry_name='Si', filename='Crystals.dat', dabax_repository=dabax_repository)

        xxraylib = xraylib.Crystal_GetCrystal('Si')

        for key in xxraylib.keys():
            tmp = xxraylib[key]

            if isinstance(tmp, list):
                for i, element in enumerate(tmp):
                    print(key, i, xdabax[key][i], xxraylib[key][i])
            else:
                print(key, xdabax[key], xxraylib[key])

    #
    # f0
    #
    if True:
        #
        # test f0 data for B3+
        #
        q = numpy.array([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9])
        f0_B3plus = numpy.array([2,1.995,1.979,1.954,1.919,1.875,1.824,1.766,1.703,1.566,1.42,1.274,1.132,0.999,0.877,0.767,0.669,0.582,0.507,0.441,0.384,0.335,0.293,0.256])

        #
        # plot
        #
        from srxraylib.plot.gol import plot

        coeff_Bdot = numpy.array([])
        plot(q, f0_B3plus,
             q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 3.0, dabax_repository=dabax_repository), q),
             q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 2.8, dabax_repository=dabax_repository), q),
             xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
             legend=["B3plus original",
                     "B3plus from f0_with_fractional_charge(5,+3)",
                     "B3plus from f0_with_fractional_charge(5,+2.8)"],
             title="")

    #
    # f0 another test
    #
    if True:
        #
        # test f0 data for B3+
        #
        q = numpy.array(
            [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
             1.7, 1.8, 1.9])
        f0_B3plus = numpy.array(
            [2, 1.995, 1.979, 1.954, 1.919, 1.875, 1.824, 1.766, 1.703, 1.566, 1.42, 1.274, 1.132, 0.999, 0.877, 0.767,
             0.669, 0.582, 0.507, 0.441, 0.384, 0.335, 0.293, 0.256])

        #
        # plot
        #
        from srxraylib.plot.gol import plot

        coeff_Bdot = numpy.array([])
        plot(q, f0_B3plus,
             q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 3.0, dabax_repository=dabax_repository), q),
             q, calculate_f0_from_f0coeff(f0_with_fractional_charge(5, 2.8, dabax_repository=dabax_repository), q),
             xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
             legend=["B3plus original",
                     "B3plus from f0_with_fractional_charge(5,+3)",
                     "B3plus from f0_with_fractional_charge(5,+2.8)"],
             title="", show=1)

    #
    # f1f2 tests
    #

    if True:
        from orangecontrib.xoppy.util.xoppy_xraylib_util import f1f2_calc

        files_f1f2 = [
            "f1f2_asf_Kissel.dat",
            "f1f2_BrennanCowan.dat",
            "f1f2_Chantler.dat",
            "f1f2_CromerLiberman.dat",
            "f1f2_EPDL97.dat",
            "f1f2_Henke.dat",
            "f1f2_Sasaki.dat",
            "f1f2_Windt.dat",
        ]
        out = []
        energy = numpy.linspace(1000,28000,500)
        F=1  #             for F in [1]: # range(12)
        tmpX = f1f2_calc("Si", energy, F=F, theta=2e-3, verbose=0)
        out.append(tmpX)

        for file_f1f2 in files_f1f2:
            tmp = f1f2_calc_dabax("Si", energy, F=F, theta=2e-3, verbose=0,
                                  filename=file_f1f2, dabax_repository=dabax_repository )
            out.append(tmp)


        legend = ['xraylib']
        for file_f1f2 in files_f1f2:
            legend.append(file_f1f2)

        plot(energy, out[0],
             energy, out[1],
             energy, out[2],
             energy, out[3],
             energy, out[4],
             energy, out[5],
             energy, out[6],
             energy, out[7],
             energy, out[8],
             legend=legend,
             xlog=1,ylog=1)

        for F in range(12):
            a_xraylib = f1f2_calc("Si", 10000, F=F, theta=2e-3, verbose=0)
            a_dabax = f1f2_calc_dabax("Si", 10000, F=F, theta=2e-3, verbose=0,
                                  filename=file_f1f2, dabax_repository=dabax_repository)
            diff = (numpy.array(a_dabax) - numpy.array(a_xraylib)) / numpy.array(a_xraylib)
            print("dabax: ", file_f1f2, a_dabax)
            print("xraylib: ", a_xraylib)
            print("diff: ", numpy.abs( diff.sum()))
            assert (numpy.abs( diff.sum()) < 0.11 )


    if True:
        #
        # misc
        #
        print("Ge, Si: ", atomic_weights_dabax(["Ge","Si"],dabax_repository=dabax_repository))
        print("70Ge: ", atomic_weights_dabax("70Ge",dabax_repository=dabax_repository))

        print(atomic_symbols_dabax()[14], atomic_names_dabax()[14])

        print("Si atomic mass", atomic_constants_dabax("Si", return_item=2, dabax_repository=dabax_repository, verbose=0))
        print("Si,Ge atomic mass", atomic_constants_dabax(["Si","Ge"], return_item=2, dabax_repository=dabax_repository, verbose=0))
        print("Si,Co atomic mass", atomic_constants_dabax(["Si", "Co"], return_label='AtomicMass', dabax_repository=dabax_repository, verbose=0))

        print("Z=27", atomic_symbols_dabax()[27])
        print("Ge Z=%d" % atomic_number_dabax("Ge"))

        print("Density Si: ", xraylib.ElementDensity(14), element_density_dabax("Si", dabax_repository=dabax_repository,verbose=0))

        # TODO: does not work for double parenthesis "Ga2(F(KI))3"
        for descriptor in ["H2O","Eu2H2.1O1.3","PO4", "Ca5(PO4)3.1F"]:
            print("\n",descriptor, dabax_repository)
            print("DABAX: ", CompoundParser(descriptor, dabax_repository=dabax_repository, verbose=0))
            print("XRAYLIB: ", xraylib.CompoundParser(descriptor))

    if True:
        #
        # cross sections
        #
        from orangecontrib.xoppy.util.xoppy_xraylib_util import cross_calc

        unit = 1
        filenames = ["CrossSec_EPDL97.dat",
                     "CrossSec_BrennanCowan.dat",
                     "CrossSec_McMaster.dat",
                     "CrossSec_NISTxaamdi.dat",
                     "CrossSec_PE_Scofield.dat",
                     "CrossSec_StormIsrael.dat",
                     "CrossSec_XCOM.dat",
                     ]

        energy = numpy.linspace(10000, 20000, 200)
        out = []
        tmpX = cross_calc("Si", energy, calculate=0, unit=unit)
        out.append(tmpX)
        for filename in filenames:
            tmp = cross_calc_dabax("Si", energy, partial='TotalCrossSection', unit=unit,
                                              filename=filename, dabax_repository=dabax_repository, verbose=0)
            out.append(tmp)

        legend = ['xraylib']
        for file in filenames:
            legend.append(file)

        plot(energy, out[0],
             energy, out[1],
             energy, out[2],
             energy, out[3],
             energy, out[4],
             energy, out[5],
             energy, out[6],
             energy, out[7],
             legend=legend,
             xlog=1,ylog=1)