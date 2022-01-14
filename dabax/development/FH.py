import numpy
from f0_util import f0_with_fractional_charge
from tmp.dabax_util import Crystal_GetCrystal, calculate_f0_from_f0coeff, Bragg_angle
import scipy.constants as codata


def Crystal_F_H_StructureFactor(crystal_id,
                                     energy_in_kev,
                                     millerH,
                                     millerK,
                                     millerL,
                                     debyeWaller,
                                     ratio_theta_thetaB=1.0):

    energy = energy_in_kev * 1e3
    wavelength = codata.h * codata.c / codata.e / energy * 1e10
    for key in crystal_id:
        print(key)
    print(crystal_id["n_atom"])
    for i in range(crystal_id["n_atom"]):
        atom_i = crystal_id["atom"][i]
        print(atom_i)
        coeffs = f0_with_fractional_charge(atom_i['Zatom'], charge=atom_i['charge'],
                                         filename="f0_InterTables.dat",
                                         dabax_repository=dabax_repository)

        angle = Bragg_angle(crystal_id, energy_in_kev,
                            millerH, millerK, millerL) * ratio_theta_thetaB
        ratio = numpy.sin(angle) / wavelength
        print(">>> ratio: ", ratio)
        f0_i = calculate_f0_from_f0coeff(coeffs, ratio)


        print(">>>>>", atom_i['Zatom'], f0_i)



if __name__ == "__main__":
    # redefine the default server at ESRF because default server name is different outside and inside ESRF
    import socket
    if socket.getfqdn().find("esrf") >= 0:
        dabax_repository = "http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"
        # dabax_repository = "/scisoft/DABAX/data"
    else:
        dabax_repository = "http://ftp.esrf.eu/pub/scisoft/DabaxFiles/"

    si_id =  Crystal_GetCrystal(entry_name='Si', filename='Crystals.dat',
                           dabax_repository=dabax_repository, verbose=True)

    F0 =  Crystal_F_H_StructureFactor(si_id,
                                    8.0,
                                    0,
                                    0,
                                    0,
                                    1.0,
                                    ratio_theta_thetaB=1.0)
    print("F0: ", F0)