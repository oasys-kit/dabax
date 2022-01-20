import numpy
from dabax.dabax_base import DabaxBase
from dabax.dabax_xoppy_decorator import DabaxXoppyDecorator

class DabaxXoppy(DabaxBase, DabaxXoppyDecorator):
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
    if socket.getfqdn().find("esrf") >= 0:
        dx = DabaxXoppy(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    else:
        dx = DabaxXoppy()

    print(dx.info())

    # f0
    if True:
        from orangecontrib.xoppy.util.xoppy_xraylib_util import f0_calc

        Si_xrl = f0_calc      (0, "Si", 0, 6, 100)
        Si_dbx = dx.f0_calc(0, "Si", 0, 6, 100)

        H2O_xrl = f0_calc      (1, "H2O", 0, 6, 100)
        H2O_dbx = dx.f0_calc(1, "H2O", 0, 6, 100)

        # TODO: not yet working with NIST
        # H2O_xrl = f0_calc      (2, "Water, Liquid", 0, 6, 100)
        # H2O_dbx = dx.f0_calc(2, "Water, Liquid", 0, 6, 100)


        from srxraylib.plot.gol import plot
        plot(Si_xrl["data"][0,:],Si_xrl["data"][1,:],
             Si_dbx["data"][0,:],Si_dbx["data"][1,:],
             H2O_xrl["data"][0, :], H2O_xrl["data"][1, :],
             H2O_dbx["data"][0, :], H2O_dbx["data"][1, :],
             linestyle=[None,'',None,''],
             marker=[None,'+',None,'+'],
             color=['r','r','b','b'],
             legend=['Si xraylib','Si dabax','H2O xraylib','H2O dabax'])





