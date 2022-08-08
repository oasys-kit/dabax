# See also the examples are in the main of
#
#  dabax_xraylib.py
#  dabax_xoppy.py
#

import numpy
from dabax.dabax_xraylib import DabaxXraylib

if __name__ == "__main__":
    dx = DabaxXraylib()

    if True:
        #
        # dabax vs xraylib tests
        #

        #
        # crystal tests
        #

        print("DABAX crystal list: \n",        dx.Crystal_GetCrystalsList())

        siD =      dx.Crystal_GetCrystal('Si')

        print("DABAX crystal si: \n",        dx.Crystal_GetCrystal('Si'))

        print("Si 111 d-spacing: DABAX: %g"% \
              (dx.Crystal_dSpacing(siD,1,1,1)))

        print("Si 111 bragg angle at 10 keV [deg]: DABAX: %g"% (
              180 / numpy.pi * dx.Bragg_angle(siD,10, 1,1,1) ))

        # F0 = dx.Crystal_F_H_StructureFactor(siD,8.0,0,0,0,1.0,ratio_theta_thetaB=1.0)
        dabax_all_F = dx.Crystal_F_0_F_H_F_H_bar_StructureFactor(siD,8.0,1,1,1,1.0,rel_angle=1.0)

        print("F0 dabax, xraylib: ",
              dx.Crystal_F_H_StructureFactor(siD,8.0,0,0,0,1.0,1.0), dabax_all_F[0])

        print("F111 dabax, xraylib: ",
              dx.Crystal_F_H_StructureFactor     (siD,8.1,1,1,1,1.0,1.0), dabax_all_F[1])

        print("F-1-1-1 dabax, xraylib: ",
              dx.Crystal_F_H_StructureFactor     (siD,8.1,-1,-1,-1,1.0,1.0), dabax_all_F[2])

        #
        # basic tools
        #
        # TODO: does not work for double parenthesis "Ga2(F(KI))3"
        for descriptor in ["H2O","Eu2H2.1O1.3","PO4", "Ca5(PO4)3.1F"]:
            print("\ncompound parsing for %s" % descriptor)
            print("DABAX: ",        dx.CompoundParser(descriptor))



        #
        # NIST compounds
        #

        nist_d = dx.GetCompoundDataNISTList()

        for i in range(len(nist_d)):
            print (dx.GetCompoundDataNISTByName(nist_d[i]))


        for i in range(len(nist_d)):
            print(nist_d[i], nist_d[i])
            a_d =      dx.GetCompoundDataNISTByIndex(i)
            print("\n\n\n", i, "\n", a_d)
            name_d = a_d["name"]

            print(name_d)
            b_d =      dx.GetCompoundDataNISTByName(name_d)
            print("\n", b_d)


        print(dx.CompoundParserCheckingNIST("H2O"))
        print(dx.CompoundParserCheckingNIST("Water, Liquid"))

        #
        # scattering factors
        #



        # loops
        energies = numpy.linspace(15,18,10)
        f1f2_d = numpy.array(dx.FiAndFii(14,energies))
        print(f1f2_d.shape)
        for i,energy in enumerate(energies):
            print("energy = %g" %energy)
            print("   Fi  dabax: ",  f1f2_d[0,i])
            print("   Fii dabax: ",  f1f2_d[1,i])


        from dabax.dabax_files import dabax_f1f2_files
        for file in dabax_f1f2_files():
            dx1 = DabaxXraylib(file_f1f2=file)
            print("\nFi  dabax (%s): %g" % (file, dx1.Fi(14, 18.0) ))
            print("Fii  dabax (%s): %g" % (file, dx1.Fii(14, 18.0) ))


        #
        # cross sections
        #

        from dabax.dabax_files import dabax_crosssec_files

        for file in dabax_crosssec_files():
            dx1 = DabaxXraylib(file_CrossSec=file)
            print(">>>>>>>>>>>>>>>  file: ", file)
            try:
                print("CSb_Total Si dabax: ",  dx1.CSb_Total(14, 18.0))
                print("CSb_Photo Si dabax: ",  dx1.CSb_Photo(14,18.0))
                print("CSb_Rayl  Si dabax: ",  dx1.CSb_Rayl (14,18.0))
                print("CSb_Compt Si dabax: ",  dx1.CSb_Compt(14,18.0))

                print("CS_Total Si dabax: ",  dx1.CS_Total(14, 18.0))
                print("CS_Photo Si dabax: ",  dx1.CS_Photo(14,18.0))
                print("CS_Rayl  Si dabax: ",  dx1.CS_Rayl (14,18.0))
                print("CS_Compt Si dabax: ",  dx1.CS_Compt(14,18.0))
            except:
                print("!!!!!!!!!!Errors with file", file)


        for file in dabax_crosssec_files():
            dx1 = DabaxXraylib(file_CrossSec=file)
            print(">>>>>>>>>>>>>>>  file: ", file)

            try:
                print("CSb_Total SiO2 dabax: ",  dx1.CSb_Total_CP("SiO2", 18.0))
                print("CSb_Photo SiO2 dabax: ",  dx1.CSb_Photo_CP("SiO2",18.0))
                print("CSb_Rayl  SiO2 dabax: ",   dx1.CSb_Rayl_CP("SiO2",18.0))
                print("CSb_Compt SiO2 dabax: ",  dx1.CSb_Compt_CP("SiO2",18.0))

                print("CS_Total SiO2 dabax: ",  dx1.CS_Total_CP("SiO2", 18.0))
                print("CS_Photo SiO2 dabax: ",  dx1.CS_Photo_CP("SiO2",18.0))
                print("CS_Rayl  SiO2 dabax: ",   dx1.CS_Rayl_CP("SiO2",18.0))
                print("CS_Compt SiO2 dabax: ",  dx1.CS_Compt_CP("SiO2",18.0))
            except:
                print("!!!!!!!!!!Errors with file", file)


        #
        # refractive index
        #
        dens = dx.element_density("Be")
        print("Refractive_Index_Re Be  dabax: ",  dx.Refractive_Index_Re ("Be",18.0, dens))
        print("Refractive_Index_Im Be  dabax: ",  dx.Refractive_Index_Im ("Be",18.0, dens))

        # loops
        energies = numpy.linspace(15,18,10)
        r1 = dx.Refractive_Index_Re("Be", energies, dens)
        for i,energy in enumerate(energies):
            print("   delta @ %g keV, dabax: %g" % (energy, 1-r1[i]))

    else:
        pass
