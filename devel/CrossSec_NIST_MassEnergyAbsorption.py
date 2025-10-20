import xraylib
import time
from scipy.constants import N_A, Avogadro

local_time = time.localtime()
print(local_time)
print(time.strftime("%Y-%m-%d %H:%M:%S", local_time))
filename = "CrossSec_NIST_MassEnergyAbsorption.dat"
f = open(filename, 'w')

header = ""
header += '#F %s' % filename
header += '\n#C  This file has been created using CrosssSec_NIST_MassEnergyAbsorption.py on ' + time.strftime("%Y-%m-%d %H:%M:%S", local_time)
header += '\n#C  Note: for Mass Attenuation coefficient and cross section see file CrossSec_NIST.dat'
header += '\n#C  This file belongs to the DABAX library. More information on'
header += '\n#C  DABAX can be found at:'
header += '\n#C  https://github.com/oasys-kit/DabaxFiles/tree/main'
header += '\n#D  ' + time.strftime("%Y-%m-%d %H:%M:%S", local_time)
header += '\n#UT Cross Sections from NIST Standard Reference Database 126'
header += '\n#UD Cross Sections from NIST Standard Reference Database 126'
header += '\n#UD'
header += '\n#UD Tables of X-Ray Mass Attenuation Coefficients and Mass '
header += '\n#UD Energy-Absorption Coefficients from 1 keV to 20 MeV for '
header += '\n#UD Elements Z = 1 to 92'
header += '\n#UD '
header += '\n#UD Reference:'
header += '\n#UD    Hubbell, J. H. and Seltzer, S. M. (1997) '
header += '\n#UD    Tables of X-Ray Mass Attenuation Coefficients and Mass '
header += '\n#UD    Energy-Absorption Coefficients (version 1.03), [Online]. '
header += '\n#UD Available: '
header += '\n#UD https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients [2025, October 20] '
header += '\n#UD National Institute of Standards and Technology, Gaithersburg, MD.'
header += '\n#UD '
header += '\n#UD Originally published as:'
header += '\n#UD X-Ray Mass Attenuation Coefficients'
header += '\n#UD NIST Standard Reference Database 126'
header += '\n#UD https://dx.doi.org/10.18434/T4D01F'
header += '\n#UD'
header += '\n#UD The  data in this file has been downloaded from: '
header += '\n#UD http://physics.nist.gov/PhysRefData/XrayMassCoef/tab3.html'
header += '\n#UD '
header += '\n#UD The data in this file are valid in the energy range 1.0 keV to 100 MeV'
header += '\n#UD for the elements Z = 1 - 100.'
header += '\n#UD'
header += '\n#UD Column description (in brackets, the corresponding label in the published tables):'
header += '\n#UD 1:PhotonEnergy[MeV]'
#header += '\n#UD 2:MassAttenuation[cm^2/g]'
header += '\n#UD 2:MassAttEnergyAbsorption[cm^2/g]'
header += '\n#UD 3:TotalCrossSection[barn/atom] (Mass Energy Absorption)'
#header += '\n#UD 5:TotalCrossSectionEnergyAborption[barn/atom]'
header += '\n#UD'
header += '\n#UD The column 1,2 are copied from the above reference. The column 3'
header += '\n#UD is been calculated from columns 2 for consistency with'
header += '\n#UD other DABAX cross section tables and related applications. The formula used is'
header += '\n#UD col3 = col2 / cf'
header += '\n#UD cf = 1e-24*AVOGADRO/AtomicMass '
header += '\n#UD For info, the keyword #UCONV contains the used numerical value 1/cf for each Z'
header += '\n#UD For the absorption edges, the energy value is duplicated (before and after the jump)'
header += '\n#UD'

f.write(header)

with open('list', 'r') as file:
    lst = file.readlines()

for i, file in enumerate(lst):
    Z = i + 1
    symbol = xraylib.AtomicNumberToSymbol(Z)
    cf = 1e-24 * Avogadro / xraylib.AtomicWeight(Z)


    file1 = file.rstrip()
    with open(file1, 'r') as file:
        text = file.readlines()

    # for j in range(len(text)): print(text[j])

    f.write('\n\n#S  %d  %s' % (i + 1, symbol))
    f.write('\n#N 5')
    f.write('\n#UCONV %.6g' % (1 / cf))
    # f.write('\n#L  PhotonEnergy[MeV]  MassAttenuation[cm^2/g]  ' + \
    #      'MassAttEnergyAbsorption[cm^2/g]  ' + \
    #      'TotalCrossSection[barn/atom]  ' + \
    #      'TotalCrossSectionEnergyAborption[barn/atom]')
    f.write('\n#L  PhotonEnergy[MeV]  ' + \
         'MassAttEnergyAbsorption[cm^2/g]  ' + \
         'TotalCrossSection[barn/atom]  ')

    for j in range(len(text)):
        line = text[j].rstrip()
        to_copy = 1
        if 'html' in line: to_copy = 0
        if 'hr>' in line: to_copy = 0
        if 'head>' in line: to_copy = 0
        if '<body' in line: to_copy = 0
        if '</body' in line: to_copy = 0
        if '<pre' in line: to_copy = 0
        if 'dl>' in line: to_copy = 0
        if '<i>' in line: to_copy = 0
        if 'Energy' in line: to_copy = 0
        if '<img' in line: to_copy = 0
        if '(MeV)' in line: to_copy = 0
        if '(MeV)' in line: to_copy = 0
        # if name in line: to_copy = 0
        for v in ['A', 'I', 'O', 'U']:
            if v in line: to_copy = 0

        line = line.replace("K",  " ")
        line = line.replace("L1", "  ")
        line = line.replace("L2", "  ")
        line = line.replace("L3", "  ")
        line = line.replace("M1", "  ")
        line = line.replace("M2", "  ")
        line = line.replace("M3", "  ")
        line = line.replace("M4", "  ")
        line = line.replace("M5", "  ")
        line = line.replace("N1", "  ")
        line = line.replace("N2", "  ")
        line = line.replace("N3", "  ")
        if not line: to_copy = 0

        if to_copy:
            tmp = line.split()
            tmp0 = float(tmp[0])
            tmp1 = float(tmp[1])
            tmp2 = float(tmp[2])
            tmp3 = tmp1 / cf
            tmp4 = tmp2 / cf

            print(line)
            # f.write("\n %.6E  %.3E  %.3E  %.6E  %.6E" % (tmp0, tmp1, tmp2, tmp3, tmp4))
            f.write("\n %.6E  %.3E  %.6E" % (tmp0, tmp2, tmp4))

f.write("\n")
f.close()
print("File %s written to disk" % filename)