"""This script going to scrape information from the CHARMM-GUI Archive for the small molecule (https://www.charmm-gui.org/?doc=archive&lib=csml).
This page contains the resname, the common name, a link to the structural view and to download the PDB for the selected lipid.
Each lipid are sorted by categories.

Usage
-----
    python scrap_charmm_gui_CMSL.py
"""

__authors__ = ("Karine DUONG", "Pierre POULAIN")
__contact__ = "pierre.poulain@u-paris.fr"


from bs4 import BeautifulSoup, Tag
import os
import pandas as pd
from pathlib import Path
import requests

from collections import Counter
import MDAnalysis as mda


def get_formula_res_name_from_pdb_file(link: str) -> tuple[str, str]:
    """Fetches a PDB file from a given link, analyzes it to extract the formula and residue name.

    Ressources
    ----------
    https://docs.mdanalysis.org/2.0.0/documentation_pages/topology/guessers.html#guessing-elements-from-atom-names
    https://docs.mdanalysis.org/2.0.0/documentation_pages/topology/guessers.html#MDAnalysis.topology.guessers.guess_atom_type

    Parameters
    ----------
        link: str
            The URL link to the PDB file.

    Returns
    -------
        tuple(str, str)
            A tuple containing the chemical formula and the name of the residue.

    Raises
    ------
        requests.exceptions.HTTPError: If an HTTP error occurs during the request.
        Exception: For other unexpected errors.
    """
    try:
        response = requests.get(link)
        response.raise_for_status()
        filename = link.split("/")[-1]
        open(filename, "wb").write(response.content)

    # In case the given link leads to nothing (error 404)
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
        return "HTTP error", None
    except Exception as err:
        print(f"An error occurred: {err}")
        return "Error", None

    # In case, the file is empty, so 0 octet
    if os.path.getsize(filename) == 0:
        return "Empty file", None

    # Read the PDB file we just download, with MDAnalylis
    molecule = mda.Universe(filename)
    atom_names = []
    res_name = set()

    # Save in a list, all the atom names and residue name from this file.
    # Unless the atom name is H.
    for atom in molecule.atoms:
        atom_name = f"{mda.topology.guessers.guess_atom_type(atom.name)}"
        # print(mda.topology.guessers.guess_atom_element(atom.name))
        if atom_name != "H":
            atom_names.append(atom_name)
        res_name.add(atom.resname)

    atom_names = Counter(atom_names)
    sorted_atom_counts = sorted(atom_names.items())
    formula = "".join(
        f"{atom}{count}" if count > 1 else atom for atom, count in sorted_atom_counts
    )

    # To delete the file we just download
    current_dir = Path.cwd()
    absolute_file_path = current_dir / filename
    relative_path = absolute_file_path.relative_to(current_dir)
    relative_path.unlink()

    print(f"File {filename} : {formula} // {res_name}")
    return formula, res_name.pop()


def parse_lipid_from_charmm_gui() -> pd.core.frame.DataFrame:
    """Parses lipid information from the CHARMM-GUI website and saves it as a CSV file.

    Returns
    -------
        pd.core.frame.DataFrame
            A DataFrame containing the lipid information with columns:
                - "Category": The category under which the lipid is listed.
                - "Alias": The alias of the lipid extracted from the download link.
                - "Name": The name of the lipid.
                - "View_Link": The link to the structural view of this lipid
                - "PDB_Link": The link to download the PDB for this lipid.
                - "Formula": The formula wwithout hydrogens
                - "Res_name_PDB": The residue name for this lipid in PDB file
    """
    url = "https://www.charmm-gui.org/?doc=archive&lib=csml"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")

    base_url = "https://www.charmm-gui.org/"
    col_names = [
        "Category",
        "Alias",
        "Name",
        "View_Link",
        "PDB_Link",
        "Formula",
        "Res_name_PDB",
    ]
    recordings = []

    rows = soup.find_all("tbody")
    for row in rows:
        category = row.get("id")
        for tr in row.find_all("tr"):
            cells = tr.find_all("td")
            if len(cells) > 1:
                alias = cells[0].text.strip()
                name = cells[1].text.strip()
                view_structure_link = f"https://www.charmm-gui.org/?doc=visualization.ngl.archive&pdb_id={alias.lower()}&arg=csml"
                download_link = (
                    base_url + cells[4].find("a")["href"]
                    if cells[4].find("a")
                    else None
                )

                formula, resname_PDB = None, None
                if download_link != None:
                    formula, resname_PDB = get_formula_res_name_from_pdb_file(
                        download_link
                    )

                recording = {
                    "Category": category,
                    "Alias": alias,
                    "Name": name,
                    "View_Link": view_structure_link,
                    "PDB_Link": download_link,
                    "Formula": formula,
                    "Res_name_PDB": resname_PDB,
                }
                recordings.append(recording)
    df = pd.DataFrame(recordings)
    df.to_csv(
        "lipid_CHARMM_GUI_CSML.csv",
        sep=",",
        index=False,
        header=True,
        columns=col_names,
        quotechar='"',
    )
    return df


if __name__ == "__main__":
    parse_lipid_from_charmm_gui()
