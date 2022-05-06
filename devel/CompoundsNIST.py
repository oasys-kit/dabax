import xraylib

f = open("CompoundsNIST.dat", 'w')

header = """#F CompoundsNIST.dat
#C  This file has been created by CompoundsNIST.py using xraylib
#C  This file belongs to the DABAX library. 
#C  DABAX files are at: http://ftp.esrf.eu/pub/scisoft/DabaxFiles/
#UT NIST Compounds and Mixtures (formula, density and composition)
#UD NIST Compounds and Mixtures (formula, density and composition)
#UD 
#UD This file contains composition and densities of a list of usual 
#UD materials. Data comes from:
#UD 3) NIST (from A-150 TISSUE-EQUIVALENT PLASTIC to WATER, LIQUID) Ref:
#UD    https://physics.nist.gov/PhysRefData/XrayMassCoef/tab4.html
#UD 
#UD Data follows the same convention as in xraylib (https://github.com/tschoonj/xraylib/)
#UD Data is arranged in a number of lines equal to the number of 
#UD elements in the mixtute sorted by Z. 
#UD The column description is: 
#UD  1: Elements (the Z value)
#UD  2: massFractions (the weight fractions)
#UD  
#UD The additional keywords are defined (as in the xraylib compound dictionary): 
#UD    #Uname
#UD    #UnElements
#UD    #Udensity
#UD 
"""

f.write(header)

list1 = xraylib.GetCompoundDataNISTList()

for i in range(len(list1)):
    a1 = xraylib.GetCompoundDataNISTByName(list1[i])
    a2 = xraylib.GetCompoundDataNISTByIndex(i)
    print("\n\n\n", list1[i], "\n", a1, "\n", a2)

    f.write("\n\n#S %i %s" % (i, list1[i]))


    f.write("\n#Uname %s"      % (a1["name"]))
    f.write("\n#UnElements %d" % (a1["nElements"]))
    f.write("\n#Udensity %g"   % (a1["density"]))

    f.write("\n#N 2")
    f.write("\n#L Elements  massFractions")

    for j in range(a1["nElements"]):
        f.write("\n%d %g" % (a1["Elements"][j], a1["massFractions"][j]))

    f.write("\n")
f.close()
print("File CompoundsNIST.dat written to disk")