"""This file lists all molecules (amino acids, ions and solvent)
used for the identification of molecular entities.

Usage
-----
    import mol_def 
"""

import MDAnalysis as mda


# https://cdn.numerade.com/ask_images/62305c89db43434fb95d4f223037af3e.jpg
BOND_LENGTH = {
    "C-C": 1.54,
    "C=C": 1.34,
    "C=O": 1.20,
    "O-H": 0.97,
    "C-H": 1.09,
    "N-H": 1.00,
    "C-S": 1.81,
    "S-H": 1.32,
    "N-C": 1.47,
    "C=N": 1.27,
    "S-S": 2.04,
    "C-O": 1.43,
}  # in Angstrom
# hydrogenBond = 2.7-3.3 A <--> 0.2-0.3 nm
# https://www.umass.edu/microbio/chime/find-ncb/help_gb.htm


# The correspondence between 3-letter code and 1-letter code
# is taken from MDAnalysis:
# https://docs.mdanalysis.org/1.0.1/_modules/MDAnalysis/lib/util.html#convert_aa_code
# See also the list of known residue names in MDAnalysis:
# https://userguide.mdanalysis.org/stable/standard_selections.html#proteins
AMINO_ACID_DICT = mda.lib.util.inverse_aa_codes
AMINO_ACID_DICT["HSP"] = "H"
AMINO_ACID_DICT["CYSP"] = "C"


NUCLEIC_ACIDS_DNA = {"DA": 'A',
                     "DT": 'T', 
                     "DC": 'C', 
                     "DG": 'G', 

                     "DA5": 'A',
                     "DA3": 'A',
                     "DT5": 'T', 
                     "DT3": 'T', 
                     "DC5": 'C', 
                     "DC3": 'C', 
                     "DG5": 'G', 
                     "DG3": 'G', 
                    }


NUCLEIC_ACIDS_RNA ={"RA": 'A', 
                    "RU": 'U', 
                    "RC": 'C', 
                    "RG": 'G', 

                    "RA5": 'A', 
                    "RA3": 'A', 
                    "RU5": 'U', 
                    "RU3": 'U', 
                    "RC5": 'C', 
                    "RC3": 'C', 
                    "RG5": 'G', 
                    "RG3": 'G', 
                    }


NUCLEIC_ACIDS = {"DA": 'A',
                "DT": 'T', 
                "DC": 'C', 
                "DG": 'G', 

                "DA5": 'A',
                "DA3": 'A',
                "DT5": 'T', 
                "DT3": 'T', 
                "DC5": 'C', 
                "DC3": 'C', 
                "DG5": 'G', 
                "DG3": 'G', 
                
                "RA": 'A', 
                "RU": 'U', 
                "RC": 'C', 
                "RG": 'G', 

                "RA5": 'A', 
                "RA3": 'A', 
                "RU5": 'U', 
                "RU3": 'U', 
                "RC5": 'C', 
                "RC3": 'C', 
                "RG5": 'G', 
                "RG3": 'G', 
                }



# The reference for each res_name and atom_names in IONS_LIST:
# https://github.com/gromacs/gromacs/tree/main/share/top
# https://mad.ibcp.fr/explore?categories=MC%3A0006
# https://www.charmm-gui.org/?doc=open_toppar&filename=stream/misc/toppar_ions_won.str 

# List of all the possibility of ions :
# https://www.unamur.be/sciences/enligne/transition/chimie/fichesderevision/revision3/listeions.htm
IONS_LIST = [
    {"name": "Aluminium ion ", "res_name": "AL3P", "atom_names": ["AL3P"]},
    {"name": "Barium ion - CHARMM", "res_name": "BAR", "atom_names": ["BAR"]},
    {"name": "Bromide ion", "res_name": "BR", "atom_names": ["BR"]},
    {"name": "Bromide ion", "res_name": "BR-", "atom_names": ["BR"]},
    {"name": "Calcium ion", "res_name": "CA", "atom_names": ["CA"]},
    {"name": "Calcium ion", "res_name": "CA2+", "atom_names": ["CA"]},
    {"name": "Calcium ion", "res_name": "CAL", "atom_names": ["CAL"]},  # CA
    {"name": "Cadmuim II cation - CHARMM", "res_name": "CD2", "atom_names": ["CD"]},
    {"name": "Cesium ion - CHARMM", "res_name": "CES", "atom_names": ["CES"]},
    {"name": "Copper ion", "res_name": "CU", "atom_names": ["CU"]},
    {"name": "Copper ion", "res_name": "CU1", "atom_names": ["CU"]},
    {"name": "Copper ion", "res_name": "CU1+", "atom_names": ["CU1+"]},
    {"name": "Copper ion", "res_name": "CU2+", "atom_names": ["CU"]},    
    {"name": "chloride ion", "res_name": "CLA", "atom_names": ["CLA"]},  # CL
    {"name": "chloride ion", "res_name": "CL", "atom_names": ["CL"]},
    {"name": "chloride ion", "res_name": "CL-", "atom_names": ["CL"]},
    {"name": "chloride ion", "res_name": "CL", "atom_names": ["CA"]},
    {"name": "Caesium ion", "res_name": "CS", "atom_names": ["CS"]},
    {"name": "Caesium ion", "res_name": "Cs+", "atom_names": ["Cs"]},
    {"name": "Fluoride ion", "res_name": "F", "atom_names": ["F"]},
    {"name": "Fluoride ion", "res_name": "F-", "atom_names": ["F"]},
    {"name": "H2PO4 ion - CHARMM", "res_name": "H2PO", "atom_names": ["P", "O1", "O2", "O3", "O4", "H1", "H2"]},
    {"name": "H3PO4 ion - CHARMM", "res_name": "H3PO", "atom_names": ["P", "O1", "O2", "O3", "O4", "H1", "H2", "H3"]},
    {"name": "HPO4 ion - CHARMM", "res_name": "HPO4", "atom_names": ["P", "O1", "O2", "O3", "O4", "H"]},
    {"name": "H2SO4 ion - CHARMM", "res_name": "H2SO", "atom_names": ["S", "O1", "O2", "O3", "O4", "H1", "H2"]},
    {"name": "HSO4 ion - CHARMM", "res_name": "HSO4", "atom_names": ["S", "O1", "O2", "O3", "O4", "H"]},
    {"name": "Hydrogen peroxide - CHARMM", "res_name": "H2O2", "atom_names": ["HP1", "OP1", "OP2", "HP2"]},
    {"name": "Iodide ion", "res_name": "I", "atom_names": ["I"]},
    {"name": "Iodide ion", "res_name": "I-", "atom_names": ["I"]},
    {"name": "big positive ion", "res_name": "IB+", "atom_names": ["IB"]},
    {"name": "Potassium ion", "res_name": "K", "atom_names": ["K"]},
    {"name": "Potassium ion", "res_name": "K+", "atom_names": ["K"]},
    {"name": "Lithium ion", "res_name": "LI", "atom_names": ["LI"]},
    {"name": "Lithium ion", "res_name": "LI+", "atom_names": ["LI"]},
    {"name": "Lithium ion", "res_name": "LIT", "atom_names": ["LIT"]},
    {"name": "Magnesium ion", "res_name": "MG", "atom_names": ["MG"]},
    {"name": "Magnesium ion", "res_name": "MG2+", "atom_names": ["MG"]},
    {"name": "Sodium ion", "res_name": "NA", "atom_names": ["NA"]},
    {"name": "Sodium ion", "res_name": "NAs", "atom_names": ["NA"]},
    {"name": "Sodium ion", "res_name": "NA+", "atom_names": ["NA"]},
    {"name": "Ammonium", "res_name": "NH4", "atom_names": ["HZ1", "NZ", "HZ2", "HZ3", "HZ4"]},
    {"name": "Hydroxide ion", "res_name": "OH", "atom_names": ["O1"]},
    {"name": "Plutonium ion", "res_name": "PU3P", "atom_names": ["PU3P"]},
    {"name": "Potassium ion", "res_name": "POT", "atom_names": ["POT"]},  # K
    {"name": "PO4 ion - CHARMM", "res_name": "PO4", "atom_names": ["P", "O1", "O2", "O3", "O4"]},
    {"name": "Rubidium ion", "res_name": "RB", "atom_names": ["RB"]},
    {"name": "Rubidium ion - CHARMM", "res_name": "RUB", "atom_names": ["RUB"]},
    {"name": "Rubidium ion", "res_name": "Rb+", "atom_names": ["Rb"]},
    {"name": "Sodium ion", "res_name": "SOD", "atom_names": ["SOD"]},  # NA
    {"name": "SO4 ion - CHARMM", "res_name": "SO4", "atom_names": ["S", "O1", "O2", "O3", "O4"]},
    {"name": "Zinc ion", "res_name": "ZN", "atom_names": ["ZN"]},
    {"name": "Zinc ion", "res_name": "ZN2", "atom_names": ["ZN"]},
    {"name": "Zinc ion", "res_name": "ZN2+", "atom_names": ["ZN"]}, 
    {"name": "Zinc ion - CHARMM", "res_name": "ZN2", "atom_names": ["ZN2"]},

    {"name": "Calcium ion - CG model with MARTINI", "res_name": "ION", "atom_names": ["CA"]},
    {"name": "chloride ion - CG model with MARTINI", "res_name": "ION", "atom_names": ["CL-"]},
    {"name": "chloride ion - CG model with MARTINI", "res_name": "ION", "atom_names": ["CL"]},
    {"name": "Sodium ion - CG model with MARTINI", "res_name": "ION", "atom_names": ["NA+"]},
    {"name": "Sodium ion - CG model with MARTINI", "res_name": "ION", "atom_names": ["NA"]},
    {"name": "Choloneion ion - CG model with MARTINI", "res_name": "ION", "atom_names": ["NC3"]}
]


# References : 
    # https://github.com/gromacs/gromacs/blob/main/share/top/oplsaa.ff/spce.itp
    # https://github.com/gromacs/gromacs/blob/main/share/top/oplsaa.ff/tip4pew.itp
    # https://github.com/gromacs/gromacs/blob/main/share/top/oplsaa.ff/tip5pe.itp
    # https://github.com/gromacs/gromacs/blob/main/share/top/amber03.ff/urea.itp 
    # https://github.com/gromacs/gromacs/blob/main/share/top/oplsaa.ff/methanol.itp
    # https://github.com/gromacs/gromacs/blob/main/share/top/oplsaa.ff/ethanol.itp
    # https://github.com/gromacs/gromacs/blob/main/share/top/oplsaa.ff/1propanol.itp
    # https://github.com/gromacs/gromacs/blob/main/share/top/gromos43a1.ff/methanol.itp
    # https://mad.ibcp.fr/explore?categories=MC%3A0008
SOLVANTS_LIST = [{"name": "water TIP3P solvant", "res_name": "TIP3", "atom_names": ["OH2"]}, 
                 {"name": "water TIP3P solvant", "res_name": "TIP3", "atom_names": ["OW", "HW1", "HW2"]}, 
                 {"name": "water TIP3P solvant", "res_name": "SOL", "atom_names": ["OH2", "H1", "H2"]},
                 {"name": "water TIP3P/spce solvant", "res_name": "SOL", "atom_names": ["OW", "HW1", "HW2"]},  # TIP3
                 {"name": "water TIP4P/tip4pew solvant", "res_name": "SOL", "atom_names": ["OW", "HW1", "HW2", "MW"]},  # TIP4
                 {"name": "water TIP4P/tip4pew solvant", "res_name": "SOL", "atom_names": ["OW1", "HW2", "HW3", "MW4"]},  # TIP4
                 {"name": "water TIP5P/tip5pe solvant", "res_name": "SOL", "atom_names": ["OW", "HW1", "HW2", "LP1", "LP2"]},  # TIP5
                 {"name": "urea solvant", "res_name": "URE", "atom_names": ["C", "O", "N1", "H11", "H12", "N2", "H21", "H22"]},
                 {"name": "organic solvant methanol/OPLS", "res_name": "MET", "atom_names": ["C", "H", "H", "H", "OA", "HO"]},
                 {"name": "organic solvant methanol", "res_name": "MET", "atom_names": ["C", "H1", "H2", "H3", "OH", "HO"]},
                 {"name": "organic solvant methanol", "res_name": "MeOH", "atom_names": ["Me1", "O2", "H3"]},
                 {"name": "organic solvant methanol", "res_name": "MOH", "atom_names": ["C", "H", "H", "H", "OA", "HO"]},
                 {"name": "organic solvant ethanol/OPLS", "res_name": "ETH", "atom_names": ["C", "H", "H", "H", "C", "H", "H", "OA", "HO"]},
                 {"name": "organic solvant ethanol", "res_name": "ETH", "atom_names": ["C", "H1", "H2", "OH", "HO", "C1", "C1", "H3", "H4", "H5"]},
                 {"name": "organic solvant propanol/OPLS", "res_name": "POL", "atom_names": ["C", "H", "H", "H", "C", "H", "H", "C", "H", "H", "OA", "HO"]},
                 {"name": "acetonitrile", "res_name": "ACN", "atom_names": ["C1", "C2", "N"]},

                 {"name": "water - in CG model with MARTINI", "res_name": "W", "atom_names": ["W"]}, 
                 {"name": "organic solvant ethanol - CG model with MARTINI", "res_name": "EOL", "atom_names": ["COH"]}, 
                 {"name": "organic solvant ether - CG model with MARTINI", "res_name": "ETH", "atom_names": ["CO"]}, 
                 {"name": "organic solvant butane - CG model with MARTINI", "res_name": "BUT", "atom_names": ["C1"]},
                 {"name": "Acetamide - CG model with MARTINI", "res_name": "ACE", "atom_names": ["NCO"]},
                 {"name": "Acetic - CG model with MARTINI", "res_name": "ACH", "atom_names": ["OOH"]},
                 {"name": "Benzene - CG model with MARTINI", "res_name": "BENZ", "atom_names": ["R1", "R2", "R3"]},
                 {"name": "Butanol - CG model with MARTINI", "res_name": "BOL", "atom_names": ["COH"]},
                 {"name": "Chlorobenzene - CG model with MARTINI", "res_name": "CB", "atom_names": ["Cl", "R2", "R3"]},
                 {"name": "Chloroform - CG model with MARTINI", "res_name": "CLF", "atom_names": ["CX"]},
                 {"name": "Chloropropane - CG model with MARTINI", "res_name": "PRX", "atom_names": ["CX"]},
                 {"name": "Cyclohexane - CG model with MARTINI", "res_name": "CHEX", "atom_names": ["R1", "R2", "R3"]},
                 {"name": "Dodecane - CG model with MARTINI", "res_name": "DEC", "atom_names": ["C", "C", "C"]},
                 {"name": "Glycerinetrioleate - CG model with MARTINI", "res_name": "TO", "atom_names": ["GLY", "ES1", "ES2", "ES3", "C1A", "D2A", "C3A", "C4A", "C1B", "D2B", "C3B", "C4B", "C1C", "D2C", "C3C", "C4C"]},
                 {"name": "Hexadecane - CG model with MARTINI", "res_name": "HD", "atom_names": ["C", "C", "C", "C"]},
                 {"name": "Methylethylsulfide - CG model with MARTINI", "res_name": "MES", "atom_names": ["CS"]},
                 {"name": "Octadecane - CG model with MARTINI", "res_name": "OD", "atom_names": ["C", "C", "C", "C", "C"]},
                 {"name": "Octane - CG model with MARTINI", "res_name": "OCT", "atom_names": ["C", "C"]},
                 {"name": "Octanol - CG model with MARTINI", "res_name": "OCO", "atom_names": ["PC", "C"]},
                 {"name": "Propane - CG model with MARTINI", "res_name": "POP", "atom_names": ["C"]},
                 {"name": "Propanol - CG model with MARTINI", "res_name": "POL", "atom_names": ["COH"]},
                 {"name": "Propanon - CG model with MARTINI", "res_name": "PON", "atom_names": ["CO"]},
                 {"name": "Propylamine - CG model with MARTINI", "res_name": "PAM", "atom_names": ["CN"]},
                 {"name": "Cis-octadecene - CG model with MARTINI", "res_name": "BUT", "atom_names": ["C", "C", "D", "C", "C"]},
                 {"name": "Butanol - in CG model", "res_name": "C4OH", "atom_names": ["COH"]}

                #  {"name": "Trans-octadecene  (in CG model)", "res_name": "ODT", "atom_names": ["C1", "C3", "C1", "C1", "C1"]},
                # https://mad.ibcp.fr/molecule/ODT?version=668651949501128350
                 ]
