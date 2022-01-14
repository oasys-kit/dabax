#
# dabax functions with the same interface as xraylib
#
import numpy
from dabax.common_tools import calculate_f0_from_f0coeff
from dabax.common_tools import parse_formula, atomic_symbols, atomic_names, atomic_number

class DabaxXoppyDecorator(object):

    def f1f2_calc(self, descriptor, energy, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True):
        raise NotImplementedError

    def f1f2_calc_mix(self, descriptor, energy, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True):
        raise NotImplementedError

    def f1f2_calc_nist(self, descriptor, energy, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True):
        raise NotImplementedError

    def cross_calc(self, descriptor, energy, calculate=0, unit=None, density=None, verbose=True):
        raise NotImplementedError

    def cross_calc_mix(self, descriptor, energy, calculate=0, unit=None, parse_or_nist=0, density=None, verbose=True):
        raise NotImplementedError

    def cross_calc_nist(self, descriptor0, energy, calculate=0, unit=None, density=None, verbose=True):
        raise NotImplementedError

    def xpower_calc(self, energies=numpy.linspace(1000.0, 50000.0, 100), source=numpy.ones(100),
                    substance=["Be"], flags=[0], dens=["?"], thick=[0.5], angle=[3.0], roughness=0.0,
                    output_file=None):
        raise NotImplementedError

    def bragg_calc(self, descriptor="Si", hh=1, kk=1, ll=1, temper=1.0,
                   emin=5000.0, emax=15000.0, estep=100.0, fileout=None):
        raise NotImplementedError

    #################################################################################
    # f0 tools
    #################################################################################
    def f0_calc(self,
        MAT_FLAG,
        DESCRIPTOR,
        GRIDSTART,
        GRIDEND,
        GRIDN,
        FILE_NAME="",
        charge=0.0,
        ):

        qscale = numpy.linspace(GRIDSTART, GRIDEND, GRIDN)

        f0 = numpy.zeros_like(qscale)

        if MAT_FLAG == 0: # element
            descriptor = DESCRIPTOR
            # for i,iqscale in enumerate(qscale):
            Z = atomic_number(descriptor)
            coeffs = self.f0_with_fractional_charge(Z, charge=charge)
            f0 = calculate_f0_from_f0coeff(coeffs, qscale)
        elif MAT_FLAG == 1: # formula
            tmp = self.compound_parser(DESCRIPTOR)
            zetas = tmp["Elements"]
            multiplicity = tmp["nAtoms"]
            for j,jz in enumerate(zetas):
                coeffs = self.f0_with_fractional_charge(jz, charge=charge)
                f0 += multiplicity[j] * calculate_f0_from_f0coeff(coeffs, qscale)
        elif MAT_FLAG == 2: # nist
            raise Exception(NotImplementedError) #TODO: implement

            Zarray = xraylib.GetCompoundDataNISTByName(DESCRIPTOR)
            zetas = Zarray["Elements"]
            fractions = numpy.array(Zarray["massFractions"])

            multiplicity = []
            for i in range(fractions.size):
                Z = zetas[i]
                symbol = atomic_symbols()[Z]
                multiplicity.append( fractions[i] / self.atomic_weights_dabax(symbol) )

            multiplicity = numpy.array(multiplicity)
            multiplicity /= multiplicity.min()

            atwt = 0.0
            for i in range(fractions.size):
                Z = zetas[i]
                atwt += multiplicity[i] * self.atomic_weights_dabax(symbol)

            print("f0_calc - nist: ")
            print("    Descriptor: ", DESCRIPTOR)
            print("    Zs: ", zetas)
            print("    n: ", multiplicity)
            print("    atomic weight: ", atwt)

            for j, jz in enumerate(zetas):
                coeffs = self.f0_with_fractional_charge(jz, charge=charge)
                f0 += multiplicity[j] * calculate_f0_from_f0coeff(coeffs, qscale)

        else:
            raise Exception("Not implemented")

        if FILE_NAME != "":
            with open(FILE_NAME, "w") as file:
                try:
                    file.write("#F %s\n"%FILE_NAME)
                    file.write("\n#S 1 xoppy f0 results\n")
                    file.write("#N 2\n")
                    file.write("#L  q=sin(theta)/lambda [A^-1]  f0 [electron units]\n")
                    for j in range(qscale.size):
                        # file.write("%19.12e  "%energy[j])
                        file.write("%19.12e  %19.12e\n"%(qscale[j],f0[j]))
                    file.close()
                    print("File written to disk: %s \n"%FILE_NAME)
                except:
                    raise Exception("f0: The data could not be dumped onto the specified file!\n")
        #
        # return
        #
        return {"application":"xoppy","name":"f0","data":numpy.vstack((qscale,f0)),"labels":["q=sin(theta)/lambda [A^-1]","f0 [electron units]"]}

