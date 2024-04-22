import MDAnalysis as mda

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


# The reference for each res_name and atom_names in IONS_LIST:
# https://github.com/gromacs/gromacs/tree/main/share/top
# List of all the possibility of ions :
# https://www.unamur.be/sciences/enligne/transition/chimie/fichesderevision/revision3/listeions.htm
IONS_LIST = [
    {"res_name": "AL3P", "atom_names": ["AL3P"]},
    {"res_name": "BR", "atom_names": ["BR"]},
    {"res_name": "BR-", "atom_names": ["BR"]},
    {"res_name": "CA", "atom_names": ["CA"]},
    {"res_name": "CA2+", "atom_names": ["CA"]},
    {"res_name": "CAL", "atom_names": ["CAL"]},  # CA
    {"res_name": "CU", "atom_names": ["CU"]},
    {"res_name": "CU1", "atom_names": ["CU"]},
    {"res_name": "CU1+", "atom_names": ["CU1+"]},
    {"res_name": "CU2+", "atom_names": ["CU"]},
    {"res_name": "CLA", "atom_names": ["CLA"]},  # CL
    {"res_name": "CL", "atom_names": ["CL"]},
    {"res_name": "CL-", "atom_names": ["CL"]},
    {"res_name": "CS", "atom_names": ["CS"]},
    {"res_name": "Cs+", "atom_names": ["Cs"]},
    {"res_name": "F", "atom_names": ["F"]},
    {"res_name": "F-", "atom_names": ["F"]},
    {"res_name": "I", "atom_names": ["I"]},
    {"res_name": "I-", "atom_names": ["I"]},
    {"res_name": "IB+", "atom_names": ["IB"]},
    {"name": "ion K", "res_name": "K", "atom_names": ["K"]},
    {"name": "ion K", "res_name": "K+", "atom_names": ["K"]},
    {"res_name": "LI", "atom_names": ["LI"]},
    {"res_name": "LI+", "atom_names": ["LI"]},
    {"res_name": "MG", "atom_names": ["MG"]},
    {"res_name": "MG2+", "atom_names": ["MG"]},
    {"res_name": "NA", "atom_names": ["NA"]},
    {"res_name": "NA+", "atom_names": ["NA"]},
    {"res_name": "OH", "atom_names": ["O1"]},
    {"name": "ion K", "res_name": "POT", "atom_names": ["POT"]},  # K
    {"res_name": "RB", "atom_names": ["RB"]},
    {"res_name": "Rb+", "atom_names": ["Rb"]},
    {"res_name": "SOD", "atom_names": ["SOD"]},  # NA
    {"res_name": "ZN", "atom_names": ["ZN"]},
    {"res_name": "ZN2+", "atom_names": ["ZN"]}
]


SOLVANTS_LIST = [{"name": "solvant water TIP3P", "res_name": "TIP3", "atom_names": ["OH2"]}, 
                 {"name": "solvant water TIP3P", "res_name": "SOL", "atom_names": ["OW", "HW1", "HW2"]},  # TIP3
                 {"name": "solvant water TIP4P", "res_name": "SOL", "atom_names": ["OW", "HW1", "HW2", "MW"]},  # TIP4
                 {"name": "solvant water TIP5P", "res_name": "SOL", "atom_names": ["OW", "HW1", "HW2", "LP1", "LP2"]}  # TIP5
                 ]