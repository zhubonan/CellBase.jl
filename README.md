## CellBase

[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://zhubonan.github.io/CellBase.jl)
[![codecov](https://codecov.io/gh/zhubonan/CellBase.jl/branch/master/graph/badge.svg?token=ACUER6ARXE)](https://codecov.io/gh/zhubonan/CellBase.jl)

Base code and routines for building/manipulating periodic crystal structures.
Hopefully useful for more complex projects with specific applications.

## Current Highlights

- `Cell{T,D}` supports auxiliary dimensions for hyperdimensional workflows.
- The current hyper-cell model uses:
  - a `3 x 3` periodic lattice
  - the first 3 position rows as periodic Cartesian coordinates
  - rows `4:end` as auxiliary non-periodic coordinates
- `get_scaled_positions(cell)` and `set_scaled_positions!(cell, scaled)` operate on the periodic `3 x N` subspace.
- XYZ I/O uses `ExtXYZ.jl` and round-trips auxiliary coordinates through `extra_dim_*` properties.
- 3D file formats such as POSCAR, RES, STRU, and `.cell` write only the first 3 Cartesian coordinates and warn when auxiliary dimensions are dropped.


## Goal

1. Ease to use - *don't make things unnecessary complex*. 
2. High performance - *when it is needed*.
3. Extensibility and flexibility.


## Acknowledgement

Elemental datasets are taken from the [SMACT](https://github.com/WMD-group/SMACT/tree/master/smact/data) package.
