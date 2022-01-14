"""

dabax: (dataBase for x-ray data)
       python module for processing remote files containing dabax

"""

import numpy
import os

from urllib.request import urlretrieve
from silx.io.specfile import SpecFile

from dabax.common_tools import calculate_f0_from_f0coeff
from dabax.common_tools import atomic_symbols, atomic_names, atomic_number
from dabax.common_tools import parse_formula

from scipy.optimize import curve_fit

class DabaxBase(object):
    def __init__(self,
                 dabax_repository="http://ftp.esrf.eu/pub/scisoft/DabaxFiles/",
                 file_f0="f0_InterTables.dat",
                 file_f1f2="f1f2_Windt.dat",
                 file_CrossSec = "CrossSec_EPDL97.dat",
                 ):

        self._dabax_repository = dabax_repository
        self._file_f0 = file_f0
        self._file_f1f2 = file_CrossSec
        self._file_CrossSec = file_f1f2

    def set_dabax_repository(self, repo):
        self._dabax_repository = repo

    def get_dabax_repository(self):
        return self._dabax_repository

    def set_file_f0(self, filename):
        self._file_f0 = filename

    def get_file_f0(self):
        return self._file_f0

    def set_file_f1f2(self, filename):
        self._file_f1f2 = filename

    def get_file_f1f2(self):
        return self._file_f1f2

    def set_file_CrossSec(self, filename):
        self._file_CrossSec = filename

    def get_file_CrossSec(self):
        return self._file_CrossSec

    def info(self):
        txt = "################  DABAX info ###########\n"
        txt += "dabax repository: %s\n" % self.get_dabax_repository()
        txt += "dabax f0 file: %s\n" % self.get_file_f0()
        txt += "dabax f1f2 file: %s\n" % self.get_file_f1f2()
        txt += "dabax CrossSec file: %s\n" % self.get_file_CrossSec()
        txt += "########################################\n"
        return txt


    def is_remote(self):
        if "http" in self.dabax_repository:
            return True
        else:
            return False


    #########################
    # common access tools
    #########################
    def get_dabax_file(self, filename, verbose=True):
        #
        # file exists in current directory
        #
        dabax_repository = self.get_dabax_repository()
        if os.path.exists(filename):
            if verbose: print("Dabax file exists in local directory: %s " % filename)
            return filename

        #
        # download remote file
        #
        if dabax_repository[0:3] == "htt" or dabax_repository[0:3] == "ftp":
            try:
                filepath, http_msg = urlretrieve(dabax_repository + filename,
                                                 filename=filename,
                                                 reporthook=None,
                                                 data=None)
                if verbose: print("Dabax file %s downloaded from %s" % (filepath, dabax_repository + filename))
                return filename
            except:
                raise Exception("Failed to download file %s from %s" % (filename, dabax_repository))

        #
        # file exists in local repository
        #
        f1 = os.path.join(dabax_repository, filename)
        if os.path.exists(f1):
            if verbose: print("Dabax file exists in local directory: %s " % f1)
            return f1

        raise Exception(FileNotFoundError)


    #########################
    # f0
    #########################

    def get_f0_coeffs_from_dabax_file(self, entry_name="Y3+", verbose=True):
        filename = self.get_file_f0()
        file1 = self.get_dabax_file(filename, verbose=verbose)
        sf = SpecFile(file1)

        flag_found = False

        for index in range(len(sf)):
            s1 = sf[index]
            name = s1.scan_header_dict["S"]

            if name.split('  ')[1] == entry_name:
                flag_found = True
                index_found = index

        if flag_found:
            return (sf[index_found].data)[:, 0]
        else:
            return []


    def f0_with_fractional_charge(self, Z, charge=0.0, verbose=True):
        filename = self.get_file_f0()

        symbol = atomic_symbols()[Z]

        if charge == 0.0:
            return self.get_f0_coeffs_from_dabax_file(entry_name=symbol)
        else:
            # retrieve all entries
            filename = self.get_dabax_file(filename, verbose=verbose)
            sf = SpecFile(filename)
            # retrieve all entries
            entries = []
            for index in range(len(sf)):
                s1 = sf[index]
                name = s1.scan_header_dict["S"]
                entries.append(name.split('  ')[1])

            # identify the entries that match the symbol
            interesting_entries = []
            charge_list = []
            index_list = []
            for i, entry in enumerate(entries):
                if entry.find(symbol) == 0:
                    if entry == symbol:
                        interesting_entries.append(entry)
                        charge_list.append(0.0)
                        index_list.append(i)
                    else:
                        entry2 = entry.replace(symbol, '')
                        try:
                            charge_item = int(entry2[::-1])  # convert to integer the reversed string
                            charge_list.append(charge_item)
                            interesting_entries.append(entry)
                            index_list.append(i)
                        except:
                            pass

            # retrieve coefficients from these interesting entries
            coefficient_list = []
            for i in index_list:
                coefficient_list.append((sf[i].data)[:, 0])

            return self.__f0_interpolate_coefficients(charge, interesting_entries, charge_list, coefficient_list,
                                                 verbose=verbose)

    def __f0_interpolate_coefficients(self, charge, interesting_entries, charge_list, coefficient_list, verbose=True):
        #
        # f0 data
        #

        nitems = len(interesting_entries)

        if nitems == 1:
            print("Warning: no interpolating of charge: only one value available for ", interesting_entries[0])
            return coefficient_list[0]

        charge_list_difference = []
        for i in range(nitems):
            charge_list_difference.append(charge_list[i] - charge)

        charge_list_difference = numpy.array(charge_list_difference)  # convert to numpy array

        if numpy.abs(charge_list_difference).min() == 0:
            idx = numpy.abs(charge_list_difference).argmin()
            if verbose: print("No interpolating needed: returning value for ", interesting_entries[idx])
            return coefficient_list[idx]

        # get the closer charge values, no matter of the sign

        sorted_indices = numpy.argsort(numpy.abs(charge_list_difference))
        sorted_index_0 = sorted_indices[0]
        sorted_index_1 = sorted_indices[1]
        delta_data = charge_list[sorted_index_1] - charge_list[sorted_index_0]
        delta_charge = charge - charge_list[sorted_index_0]
        delta = delta_charge / delta_data
        if verbose: print("Interpolating charge %g = %s + %g (%s - %s)" % (charge,
                                                                           interesting_entries[sorted_index_0],
                                                                           delta,
                                                                           interesting_entries[sorted_index_1],
                                                                           interesting_entries[sorted_index_0]))

        # interpolate to get the f0 for the wanted charge

        q = numpy.linspace(0.0, 2.0, 100)
        f0_0 = calculate_f0_from_f0coeff(coefficient_list[sorted_index_0], q)
        f0_1 = calculate_f0_from_f0coeff(coefficient_list[sorted_index_1], q)
        f0 = f0_0 + delta * (f0_1 - f0_0)

        #
        # fit
        #
        try:
            popt, pcov = curve_fit(self.__f0func, q, f0, p0=coefficient_list[sorted_index_0], maxfev=20000)
            if verbose: print("fitted: ", popt)

            return popt
        except:
            if verbose: print("Error: failed to fit coefficients for fractional charge. Returning the ones of ",
                              interesting_entries[sorted_index_0])
            return coefficient_list[sorted_index_0]

    def __f0func(self, q, a1, a2, a3, a4, a5, a6, a7, a8, a9):
        return calculate_f0_from_f0coeff([a1, a2, a3, a4, a5, a6, a7, a8, a9], q)

    ######################
    # miscellaneous
    ######################

    def atomic_weights(self, descriptor,
                             filename="AtomicWeights.dat",
                             verbose=True,
                             ):
        """
        ;       Returns atomic weights from DABAX.
        ;
        ; INPUTS:
        ;       id: an identifier string (i.e. 'Si', '70Ge)
        ;
        ;       If descriptor is the symbol (e.g., Ge),
        ;         the averaged atomic mass is returned.
        ;       If descriptor contains the isotope (number of nucleons) (e.g., 70Ge),
        ;         the atomic mass for the isotope is returned.
        ;
        ;       filename = the DABAX  inout file (default AtomicWeights.dat)

        """

        if isinstance(descriptor, str):
            descriptor = [descriptor]
            descriptor_is_string = 1
        else:  # is list
            descriptor_is_string = 0

        # access spec file
        file1 = self.get_dabax_file(filename, verbose=verbose)
        sf = SpecFile(file1)

        out = []

        for idescriptor in descriptor:
            flag_found = False
            index_found = []
            for index in range(len(sf)):
                s1 = sf[index]
                name = s1.scan_header_dict["S"]
                line = " ".join(name.split())
                scan_name = line.split(' ')[1]
                if scan_name[-len(idescriptor):] == idescriptor:
                    flag_found = True
                    index_found.append(index)

            if not flag_found:
                raise (Exception("Entry name not found: %s" % idescriptor))

            data = sf[index_found[0]].data

            if idescriptor[0].isdigit():
                out.append(data[0, 0])
            else:
                out.append(data[2, 0])

        if descriptor_is_string:
            return out[0]
        else:
            return out


    def atomic_constants(self, descriptor,
                         filename="AtomicConstants.dat",
                         verbose=True,
                         return_item=0,
                         return_label=None,
                         ):
        """
        ;	Returns atomic constants from DABAX.
        ;
        ; CALLING SEQUENCE:
        ;	out = atomic_constants(id,file,return=return)
        ; INPUTS:
        ;	id: an identifier (or an array of identifiers) to be found in the
        ;	scan title (i.e. 'Si')
        ;
        ; KEYWORDS:
        ;	File = the DABAX  inout file (default: AtomicConstants.dat)
        ;	return_label and return_item  define the variable to be returned.
        ;   If return_name is not None, it has priority over retirn_index
        ;		number of the column in the DABAX file, or a text
        ;		identifier (case insensitive) listed below:
        ;		return_label='AtomicRadius'	             or return_item=0
        ;		return_label='CovalentRadius'	         or return_item=1
        ;		return_label='AtomicMass'	             or return_item=2
        ;		return_label='BoilingPoint'	             or return_item=3
        ;		return_label='MeltingPoint'	             or return_item=4
        ;		return_label='Density'	                 or return_item=5
        ;		return_label='AtomicVolume'	             or return_item=6
        ;		return_label='CoherentScatteringLength'	 or return_item=7
        ;		return_label='IncoherentX-section'	     or return_item=8
        ;		return_label='Absorption@1.8A'	         or return_item=9
        ;		return_label='DebyeTemperature'          or return_item=10
        ;		return_label='ThermalConductivity'       or return_item=11
        ;
        ; OUTPUT:
        ;	out: the value of the selected parameter
        ;
        ; EXAMPLES:
        ;	print(atomic_constants('Si',return='AtomicMass'))
        ;	    28.085500
        ;	print(atomic_constants(14,return='AtomicMass'))
        ;           28.085500
        ;	print(atomic_constants([14,27],return='AtomicMass'))
        ;	    28.085500       58.933200
        ;
        ;-

        """

        if isinstance(descriptor, str):
            descriptor = [descriptor]
            descriptor_is_string = 1
        else:  # is list
            descriptor_is_string = 0

        return_index = -1
        if return_label is None:
            return_index = return_item
        else:
            if return_label == 'AtomicRadius'	: return_index = 0
            if return_label == 'CovalentRadius'	            : return_index = 1
            if return_label == 'AtomicMass'	                : return_index = 2
            if return_label == 'BoilingPoint'	            : return_index = 3
            if return_label == 'MeltingPoint'	            : return_index = 4
            if return_label == 'Density'	                : return_index = 5
            if return_label == 'AtomicVolume'	            : return_index = 6
            if return_label == 'CoherentScatteringLength'	: return_index = 7
            if return_label == 'IncoherentX-section'	    : return_index = 8
            if return_label == 'Absorption@1.8A'	        : return_index = 9
            if return_label == 'DebyeTemperature'           : return_index = 10
            if return_label == 'ThermalConductivity'        : return_index = 11


        if return_index == -1: raise Exception("Bad item index")
        # access spec file

        file1 = self.get_dabax_file(filename, verbose=verbose)

        sf = SpecFile(file1)

        out = []

        for idescriptor in descriptor:

            for index in range(len(sf)):
                flag_found = False
                s1 = sf[index]
                name = s1.scan_header_dict["S"]
                scan_name = name.split('  ')[1]
                if scan_name == idescriptor:
                    flag_found = True
                    index_found = index
                    break

            if flag_found:
                out.append((sf[index_found].data)[return_index, 0])
            else:
                raise Exception("Data not found for %s " % idescriptor)

        if descriptor_is_string:
            return out[0]
        else:
            return out


    def element_density(self, descriptor,
                        filename="AtomicConstants.dat", verbose=True, ):

        return self.atomic_constants(descriptor, filename=filename, return_label="Density",
                                     verbose=verbose)


    def compound_parser(self, descriptor, verbose=True):

        zetas, fatomic = parse_formula(formula=descriptor, verbose=verbose)

        elements = []
        atomic_weight = []
        massFractions = []

        for i ,z in enumerate(zetas):
            symbol = atomic_symbols()[z]
            atw = self.atomic_weights(symbol, verbose=verbose)
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

if __name__ == '__main__':
    import socket
    if socket.getfqdn().find("esrf") >= 0:
        dx = DabaxBase(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    else:
        dx = DabaxBase()

    print(dx.info())


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

        plot(q, f0_B3plus,
             q, calculate_f0_from_f0coeff(dx.f0_with_fractional_charge(5, 3.0), q),
             q, calculate_f0_from_f0coeff(dx.f0_with_fractional_charge(5, 2.8), q),
             xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
             legend=["B3plus original",
                     "B3plus from f0_with_fractional_charge(5,+3)",
                     "B3plus from f0_with_fractional_charge(5,+2.8)"],
             marker=['+', None, None],
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

        plot(q, f0_B3plus,
             q, calculate_f0_from_f0coeff(dx.f0_with_fractional_charge(5, 3.0), q),
             q, calculate_f0_from_f0coeff(dx.f0_with_fractional_charge(5, 2.8), q),
             xtitle=r"q (sin $\theta$ / $\lambda$)", ytitle="f0 [electron units]",
             legend=["B3plus original",
                     "B3plus from f0_with_fractional_charge(5,+3)",
                     "B3plus from f0_with_fractional_charge(5,+2.8)"],
             marker=['+',None,None],
             title="", show=1)


    if True:
        #
        # misc
        #
        print("Ge, Si: ", dx.atomic_weights(["Ge","Si"]))
        print("70Ge: ", dx.atomic_weights("70Ge"))

        print(atomic_symbols()[14], atomic_names()[14])

        print("Si atomic mass", dx.atomic_constants("Si", return_item=2, verbose=0))
        print("Si,Ge atomic mass", dx.atomic_constants(["Si", "Ge"], return_item=2, verbose=0))
        print("Si,Co atomic mass", dx.atomic_constants(["Si", "Co"], return_label='AtomicMass', verbose=0))

        print("Z=27", atomic_symbols()[27])
        print("Ge Z=%d" % atomic_number("Ge"))

        print("Density Si: ", dx.element_density("Si", verbose=0))







