"""
Utilities for listing and accessing DABAX data files.
"""

def dabax_f0_files():
    """Return the list of available DABAX f0 (form factor) data files."""
    return ["f0_InterTables.dat", "f0_WaasKirf.dat", "f0_xop.dat"]

def dabax_f1f2_files():
    """Return the list of available DABAX f1/f2 (anomalous scattering factor) data files."""
    return ["f1f2_Windt.dat",
            "f1f2_Sasaki.dat",
            "f1f2_Henke.dat",
            "f1f2_EPDL97.dat",
            "f1f2_CromerLiberman.dat",
            "f1f2_Chantler.dat",
            "f1f2_BrennanCowan.dat",
            "f1f2_asf_Kissel.dat",
            ]

def dabax_crosssec_files():
    """Return the list of available DABAX cross-section data files."""
    return ["CrossSec_EPDL97.dat",
            "CrossSec_XCOM.dat",
            "CrossSec_BrennanCowan.dat",
            "CrossSec_StormIsrael.dat",
            "CrossSec_McMaster.dat",
            "CrossSec_NIST.dat",                      # does not contain partial cross sections (compton, rayl, photo)
            "CrossSec_NIST_MassEnergyAbsorption.dat", # does not contain partial cross sections (compton, rayl, photo)
            # "CrossSec_PE_Scofield.dat",   # does not contain partial cross sections (compton, rayl, photo)
            ]
