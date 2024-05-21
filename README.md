# GROdecoder

GROdecoder is a Python module designed to extract and identify molecules from a .gro file.

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

## Usage

Run GROdecoder on a test file:
```bash
python grodecoder.py --input data/examples/barstar.gro
```

Draw graphs corresponding to the molecules found in the structure file (by default at false):
```bash
python grodecoder.py --input data/examples/barstar.gro --drawgraph
```

Check that all residues exist only one time (by default at false):
```bash
python grodecoder.py --input data/examples/barstar.gro --checkoverlapping
```

Add edges and degre in the fingerprint (by default at false)"
```bash
python grodecoder.py --input data/examples/barstar.gro --checkconnectivity
```

Choose the method to calculate the atom pairs. If we know the resolution of the system is all atom choose 'aa' or if we don't know choose 'auto' (by default at aa): 
```bash
python grodecoder.py --input data/examples/barstar.gro --bondthreshold
```




