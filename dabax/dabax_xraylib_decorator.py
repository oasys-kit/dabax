#
# dabax functions with the same interface as xraylib
#
import numpy
from dabax.dabax_base import DabaxBase
from orangecontrib.xoppy.util.xoppy_xraylib_util import bragg_metrictensor
import scipy.constants as codata
from silx.io.specfile import SpecFile
from dabax.common_tools import parse_formula, atomic_symbols_dabax, atomic_names_dabax, atomic_number_dabax

class DabaxXraylibDecorator(object):

    #########################
    # crystal
    #########################
    def Crystal_GetCrystal(self, entry_name='YB66', filename='Crystals.dat',
                           verbose=True):
        """
        parse a complex crystal structure file into a dictionary (like xraylib.Crystal_GetCrystal('Si'))

        it has an additional fiels for each atom: the charge

        return a dictionary containing crystal infomation
        """
        dabax_repository = self.get_dabax_repository()

        file1 = self.get_dabax_file(filename, verbose=verbose)

        sf = SpecFile(file1)

        flag_found = False

        for index in range(len(sf)):
            s1 = sf[index]
            name = s1.scan_header_dict["S"]
            if name.split(' ')[1] == entry_name:
                flag_found = True
                index_found = index

        if not flag_found:
            raise (Exception("Entry name not found: %s" % entry_name))

        cryst = {'name': entry_name}  # returned dictionary like that one created by xraylib.Crystal_GetCrystal(descriptor)

        cell_parameters = sf[index_found].scan_header_dict["UCELL"]
        cell_parameters = ' '.join(cell_parameters.split())  # remove multiple blanks

        a = cell_parameters.split(' ')
        cryst['a'] = float(a[0])
        cryst['b'] = float(a[1])
        cryst['c'] = float(a[2])
        cryst['alpha'] = float(a[3])
        cryst['beta'] = float(a[4])
        cryst['gamma'] = float(a[5])

        volume = bragg_metrictensor(float(a[0]), float(a[1]), float(a[2]),
                                    float(a[3]), float(a[4]), float(a[5]), RETURN_VOLUME=1)

        cryst['volume'] = volume

        cell_data = numpy.array(sf[index_found].data)

        cryst['n_atom'] = cell_data.shape[1]
        atom = []

        for i in range(cell_data.shape[1]):
            if cell_data.shape[0] == 5:  # standard 5 columns
                # not here, this info is not in the dabax file
                # s = symbol_to_from_atomic_number(int(cell_data[0,i]))
                atom.append({
                    'Zatom': int(cell_data[0, i]),
                    'fraction': cell_data[1, i],
                    'x': cell_data[2, i],
                    'y': cell_data[3, i],
                    'z': cell_data[4, i],
                    'charge': 0.0, })
            else:  # 6 columns (charge)
                # 'AtomicName' required to compatible my current code
                # s = symbol_to_from_atomic_number(int(cell_data[0,i]))
                # if cell_data[5, i] != 0:  #charged
                #     s = s + f'%+.6g'%cell_data[5, i]
                atom.append({
                    # 'AtomicName': s,
                    'Zatom': int(cell_data[0, i]),
                    'fraction': cell_data[1, i],
                    'x': cell_data[2, i],
                    'y': cell_data[3, i],
                    'z': cell_data[4, i],
                    'charge': cell_data[5, i], })

        cryst['atom'] = atom
        cryst['cpointer'] = None

        ANISO_KEY = "UANISO_COFF"  # prefix for a line with anisotropic coefficients
        d = sf[index_found].scan_header_dict
        AnisoItem = {'Name': '       ',
                     'start': 0,
                     'end': 0,
                     'beta11': 0.0,
                     'beta22': 0.0,
                     'beta33': 0.0,
                     'beta12': 0.0,
                     'beta13': 0.0,
                     'beta23': 0.0}

        a = [(x, d[x].split()) for x in d if x[:len(ANISO_KEY)] == ANISO_KEY]
        if len(a) > 0:  # found Anisotropic coefficients in the header, process it
            a = sorted(a, key=lambda x: int(x[1][0]),
                       reverse=False)  # sort 'Start' ascendant, avoid order changed by the SpecFile
            n = 0
            Aniso = []
            for x in a:  # tuple('UANISO_COFF_B1',[1 96 0.00038 0.00044 0.00039 0 0 0])
                AnisoItem['Name'] = x[0][len(ANISO_KEY) + 1:]  # get site atom name starting from 13th character 'B1', etc
                AnisoItem['start'] = int(x[1][0])
                AnisoItem['end'] = int(x[1][1])
                AnisoItem['beta11'] = float(x[1][2])
                AnisoItem['beta22'] = float(x[1][3])
                AnisoItem['beta33'] = float(x[1][4])
                AnisoItem['beta12'] = float(x[1][5])
                AnisoItem['beta13'] = float(x[1][6])
                AnisoItem['beta23'] = float(x[1][7])
                Aniso.append(AnisoItem.copy())
                n = n + 1
            cryst['Aniso'] = Aniso  # if having key 'Ansio' when there is anisotropic data,otherwise no
            cryst['n_aniso'] = n
        else:  # create a dummy Aniso to compatible with xraylib
            cryst['Aniso'] = [AnisoItem.copy()]
            cryst['n_aniso'] = 1

        return cryst


    def Crystal_GetCrystalsList(self, verbose=True):
        """
        get crystal names from crystals.dat
        """
        file1 = self.get_dabax_file('Crystals.dat', verbose=verbose)
        sf = SpecFile(file1)
        crystals = []
        for index in range(len(sf)):
            s1 = sf[index]
            name = s1.scan_header_dict["S"]
            crystals.append(name.split(' ')[1])

        return crystals

    def Crystal_dSpacing(self, cryst, h, k, l):
        return bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'],
                                  cryst['alpha'], cryst['beta'], cryst['gamma'],
                                  HKL=[h, k, l])

    def Bragg_angle(self, cryst, E_keV, h, k, l):
        dspacing = self.Crystal_dSpacing(cryst, h, k, l)  # in A
        wavelength = codata.h * codata.c / codata.e / (E_keV * 1e3) * 1e10  # in A
        return numpy.arcsin(wavelength / 2 / dspacing)
    #
    #
    #   TODO:
    #          F_0 = xraylib.Crystal_F_H_StructureFactor(_crystal, E_keV, h, k, l, _debyeWaller, 1.0)
    #
    #          F_H = xraylib.Crystal_F_H_StructureFactor(_crystal, E_keV, h, k, l, _debyeWaller, 1.0)
    #

    def CompoundParser(self, descriptor, verbose=True):

        zetas, fatomic = parse_formula(formula=descriptor, verbose=verbose)

        elements = []
        atomic_weight = []
        massFractions = []

        for i ,z in enumerate(zetas):
            symbol = atomic_symbols_dabax()[z]
            atw = self.atomic_weights_dabax(symbol, verbose=verbose)
            elements.append(z)
            atomic_weight.append(atw)
            massFractions.append(fatomic[i ] *atw)

        mweight = 0.0
        for i in range(len(fatomic)):
            mweight += atomic_weight[i] * fatomic[i]

        for i in range(len(massFractions)):
            massFractions[i] /= mweight

        new_dict = {
            "nElements": len(elements),
            "nAtomsAll": float(numpy.array(fatomic).sum()),
            "Elements" :zetas,
            "massFractions": massFractions,
            "nAtoms" :fatomic,
            "molarMass": mweight,
        }

        return new_dict

