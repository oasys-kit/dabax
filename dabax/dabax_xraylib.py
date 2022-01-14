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
    import socket
    import xraylib

    if socket.getfqdn().find("esrf") >= 0:
        dabax_repository = "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"
    else:
        dabax_repository = "http://ftp.esrf.eu/pub/scisoft/DabaxFiles/"

    dx = DabaxXraylib(dabax_repository=dabax_repository)
    print(dx.info())

    #
    # crystal tests
    #
    if True:
        print(dx.get_dabax_file("Crystals.dat", verbose=0))

        print(dx.get_f0_coeffs_from_dabax_file("Y3+"))

        print(dx.Crystal_GetCrystalsList())

        yb = dx.Crystal_GetCrystal('YB66', filename='Crystals.dat')

        si = dx.Crystal_GetCrystal("Si")
        print("Si 111 d-spacing: ", dx.Crystal_dSpacing(si,1,1,1))
        print("Si 111 bragg angle at 10 keV [deg]: ", 180 / numpy.pi * dx.Bragg_angle(si,10, 1,1,1))

    #
    # crystal vs xraylib tests
    #

    # TODO: does not work for double parenthesis "Ga2(F(KI))3"
    for descriptor in ["H2O","Eu2H2.1O1.3","PO4", "Ca5(PO4)3.1F"]:
        print("\n",descriptor)
        print("DABAX: ", dx.CompoundParser(descriptor, verbose=0))
        print("XRAYLIB: ", xraylib.CompoundParser(descriptor))