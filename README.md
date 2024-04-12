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
python grodecoder.py --gro data/examples/barstar.gro
```

Draw graphs corresponding to the molecules found in the .gro file:

```bash
python grodecoder.py --gro data/examples/barstar.gro --drawgraph
```







