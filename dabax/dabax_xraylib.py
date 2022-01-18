import numpy
from dabax.dabax_base import DabaxBase
from dabax.dabax_xraylib_decorator import DabaxXraylibDecorator

class DabaxXraylib(DabaxBase, DabaxXraylibDecorator):
    def __init__(self,
                 dabax_repository="http://ftp.esrf.eu/pub/scisoft/DabaxFiles/",
                 file_f0="f0_InterTables.dat",
                 file_f1f2="f1f2_Windt.dat",
                 file_CrossSec = "CrossSec_EPDL97.dat",
                 ):

        DabaxBase.__init__(self,
                           dabax_repository=dabax_repository,
                           file_f0=file_f0,
                           file_f1f2=file_f1f2,
                           file_CrossSec=file_CrossSec)

if __name__ == "__main__":
    import xraylib

    import socket
    if socket.getfqdn().find("esrf") >= 0:
        dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    else:
        dx = DabaxXraylib()

    print(dx.info())

    #
    # crystal tests
    #

    print("DABAX crystal list: \n",        dx.Crystal_GetCrystalsList())
    print("XRAYLIB crystal list: \n", xraylib.Crystal_GetCrystalsList())

    siD = dx.Crystal_GetCrystal('Si')
    siX = xraylib.Crystal_GetCrystal('Si')

    print("DABAX crystal si: \n",        dx.Crystal_GetCrystal('Si'))
    print("XRAYLIB crystal si: \n", xraylib.Crystal_GetCrystal('Si'))

    print("Si 111 d-spacing: DABAX: %g, XRAYLIB: %g "% \
          (dx.Crystal_dSpacing(siD,1,1,1),dx.Crystal_dSpacing(siX,1,1,1)))

    print("Si 111 bragg angle at 10 keV [deg]: DABAX: %g, XRAYLIB: %g "% (\
          180 / numpy.pi * dx.Bragg_angle(siD,10, 1,1,1), \
          180 / numpy.pi * dx.Bragg_angle(siX, 10, 1, 1, 1)))

    #
    # dabax vs xraylib tests
    #

    # TODO: does not work for double parenthesis "Ga2(F(KI))3"
    for descriptor in ["H2O","Eu2H2.1O1.3","PO4", "Ca5(PO4)3.1F"]:
        print("\ncompound parsing for %s" % descriptor)
        print("DABAX: ", dx.CompoundParser(descriptor))
        print("XRAYLIB: ", xraylib.CompoundParser(descriptor))

    print("Si is Z= %d (DABAX)  %d (XRAYLIB)" % (dx.SymbolToAtomicNumber("Si"),xraylib.SymbolToAtomicNumber("Si")))
    print("Z=23 is %s (DABAX)  %s (XRAYLIB)" % (dx.AtomicNumberToSymbol(23),xraylib.AtomicNumberToSymbol(23)))
    print("Density Z=30 %g (DABAX)  %g (XRAYLIB)" % (dx.ElementDensity(30),xraylib.ElementDensity(30)))
    print("AtWeight Z=30 %g (DABAX)  %g (XRAYLIB)" % (dx.AtomicWeight(30),xraylib.AtomicWeight(30)))


    # crystal structure factors

    print("Fi  dabax,xraylib: ",  dx.Fi (14,18.0), xraylib.Fi (14,18.0))
    print("Fii dabax,xraylib: ", dx.Fii(14,18.0), xraylib.Fii(14,18.0))

    # F0 = dx.Crystal_F_H_StructureFactor(siD,8.0,0,0,0,1.0,ratio_theta_thetaB=1.0)

    dabax_all_F = dx.Crystal_F_0_F_H_F_H_bar_StructureFactor(siD,8.0,1,1,1,1.0,rel_angle=1.0)

    print("F0 dabax, xraylib: ",
          dx.Crystal_F_H_StructureFactor(siD,8.0,0,0,0,1.0,1.0), dabax_all_F[0],
          xraylib.Crystal_F_H_StructureFactor(siX,8.0,0,0,0,1.0,1.0))

    print("F111 dabax, xraylib: ",
          dx.Crystal_F_H_StructureFactor     (siD,8.1,1,1,1,1.0,1.0), dabax_all_F[1],
          xraylib.Crystal_F_H_StructureFactor(siX,8.1,1,1,1,1.0,1.0))

    print("F-1-1-1 dabax, xraylib: ",
          dx.Crystal_F_H_StructureFactor     (siD,8.1,-1,-1,-1,1.0,1.0), dabax_all_F[2],
          xraylib.Crystal_F_H_StructureFactor(siX,8.1,-1,-1,-1,1.0,1.0))

    print("FF_rayl dabax, xraylib: ",
          dx.FF_Rayl(17, 2.2),xraylib.FF_Rayl(17, 2.2) )