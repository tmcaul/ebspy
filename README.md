# EBSPy

This is a simple package for handling EBSPs in python. 

Currently contains:
- bImportEBSPs, load an EBSD dataset from a Bruker Esprit exported HDF5 file.
- EBSP class, a stack of EBSD data with a few useful methods. Can feed a column vector with metadata or a rectangular array to construct dataset. Add-able, edit-able, set-able using magic methods.
- EBSP.bgcorr() class method, adapted from AstroEBSD https://github.com/benjaminbritton/AstroEBSD MATLAB BGCor function.
- Varimax and generalised orthogonal rotation implementations

Background correction settings:
- Square crop (default true)
- Radial masking (default false)
- Split BG chip (default false)
- Gaussian filter EBSP (default true)
- Gaussian sigma (default 4)
- Line-fix for Bruker EBSPs (default true)
- Mean centre and stdev normalise patterns (default true)
- NMF normalise ie. make positive (default false)
- resize patterns (default false)
- EBSP screensize (default [200,200])


To install this package using anaconda:
- Download from github https://github.com/tmcaul/EBSPy
- Navigate to the package's folder using an anaconda terminal
- Run the command 'pip install .'
- Associated dependency packages should download automatically