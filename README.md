# GROdecoder

GroDecoder extracts and identifies the molecular components of a structure file (PDB or GRO) issued from a molecular dynamics simulation. 

## Installation

Clone the project:

```bash
git clone https://github.com/pierrepo/grodecoder.git
cd grodecoder
```

Create and activate a conda environment:

```bash
conda env create -f environment.yml
conda activate grodecoder-env
```

## Usage as command line tool

Run GROdecoder on a test file:
```bash
python grodecoder.py --input data/examples/barstar.gro
```

Add edges and degre in the fingerprint (by default at false)"
```bash
python grodecoder.py --input data/examples/barstar.gro --checkconnectivity
```

Choose the method to calculate the atom pairs. If we know the resolution of the system is coarse-grain enter a threshold (a positiv float number) or we don't know so choose 'auto' (by default at 'auto'): 
```bash
python grodecoder.py --input data/examples/barstar.gro --bondthreshold [auto or a threshold]
```

Add PDB id, their putative name and the organism name in the JSON file for each protein sequence (by default at false):
```bash
python grodecoder.py --input data/examples/barstar.gro --querypdb
```

## Run the web app

Run the Streamlit web app:

```bash
streamlit run streamlit_app.py
```

then open your web browser at <http://localhost:8501>

or with the URL:
```bash
https://grodecoder.streamlit.app/
```

## Update the database

Run the script for coarse-grain model molecule data: 
```bash
python script/scrap_MAD.py
```

And run the script for all-atom model molecule data: 
```bash
python script/scrap_charmm_gui_CSML.py
```

These scripts analyzes the database MAD (https://mad.ibcp.fr/explore) and CHARMM-GUI CSML (https://www.charmm-gui.org/?doc=archive&lib=lipid). It download the data - if it's not already exist. It extracts information that can be used to identify lipids in the GRO and PDB files and save it into a CSV file.