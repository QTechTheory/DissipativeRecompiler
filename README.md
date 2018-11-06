# DissipativeRecompiler

This repository contains C code for 3 executables.

- `simulator.c` which emulates real-time simulation using the [Li method](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.7.021050), Trotter method, and direct Hamiltonian exponentiation.
- `recompiler.c` which uses [imaginary-time simulation](https://arxiv.org/abs/1804.03023) to recompile one quantum circuit into another.
- `eliminator.c` which uses constrained imaginary time simulation to compactify a quantum circuit.

These can be compiled using [GNU Make](https://www.gnu.org/software/make/) with the commands
```bash
make -f makefile_simulator
make -f makefile_recompiler
make -f makefile_eliminator
```
though this requires installing the GSL as outlined [here](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903).

## Code

The repository includes a copy of [QuEST](https://quest.qtechtheory.org/) for simulating quantum circuits, in addition to

- `circuitloader.c`

> to be continued

## Usage
