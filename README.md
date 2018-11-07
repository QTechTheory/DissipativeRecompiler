DissipativeRecompiler
=====================

> **This README is a work in progress!**

**Table of Contents**
  - [Overview](#overview)
  - [Simulator](#simulator)
  - [Recompiler](#recompiler)
  - [Eliminator](#eliminator)
  - [Demonstration](#demonstration)

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

------------------------------------------------------------------------------------------------

## Overview

The repository includes a copy of [QuEST](https://quest.qtechtheory.org/) for simulating quantum circuits, in addition to

- `circuitloader.c` which parses circuit files and applies circuits to states
- `hamiltonianloader.c` which performs energy calculations
- `paramevolver.c` which evolves circuit parameters according to Li's variational method
- `trotterevolver.c` which produces states via Trotter's method
- `trueevoler.c` which produces states via direct unitary time evolution
- `utilities.c` containing simple array and file parsing utilities
- `mmaformatter.c` used to write arrays to files readable by Mathematica

------------------------------------------------------------------------------------------------


## Simulator

`simulator.c` emulates realtime simulation using Li's variational algorithm, and compares its performance to simulation via the Trotter method and direct unitary time evolution. As input, we take a Hamiltonian which is a weighted sum of Pauli operators, some parameterised circuit, and some starting wavefunction.
Li's algorithm then finds parameters which produce wavefunctions which approximate future evolution states of the starting wavefunction under the Hamiltonian.

### How it works

We take some initial state and an ansatz with initially (very near) zero parameters which effects the identity. Every future assignment of the parameters is calculated using Li's variational equations. These equations are based off expected values of the (derivatives of the) wavefunction which the ansatz produces using the current parameters, and the system Hamiltonian. 

On quantum hardware, these expected values are averages of measurements on states produced by variations of the ansatz circuit. Our simulation code instead computes these quantities directly by approximating the wavefunction derivatives using finite difference formulae. These expected values populate a family of linear equations, which we solve using GSL's Tikhonov regularisation functions. This occurs in `paramevolver.c`. The solution informs our next choice of parameters, which produce the next time-dependent state of the simulated system. We monitor the energy of the simulated state (using `hamiltonianloader.c`) and its fidelity against the "true" simulation state (produced by `trueevolver.c`) and that produced by the Trotter method (`trotterevolver.c`).

Note that `trueevolver.c` still approximates the evolution, but at a substantially greater accuracy than the Li and Trotter approaches. It achieves this by re-expressing the Pauli Hamiltonian in the computational basis, and numerically exponentiating this matrix (using GSL) to build the unitary time evolution operator for some timestep. This is matrix multiplied with a wavefunction's statevector to evolve it forward in time by the duration timestep.

Meanwhile, `trotterevolver.c` accepts a fixed-length Trotter circuit and chooses the gate parameters (by Trotterisation) in order to immediately produce a future state. It can optionally randomise the gate order or reverse every second Trotter cycle in order to improve its accuracy.

### Usage

After compiling with
```bash
make -f makefile_simulator
```
the simulator is run by suppling command line arguments
```bash
./simulator ansatz in_params out_params out_truewavef out_data [in_truewavef]
````
where

| argument       | explanation |
|----------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| `ansatz`         | filename of the input ansatz circuit                                                                                                           |
| `in_params`      | filename of the input initial ansatz parameters                                                                                                |
| `out_params`     | filename to which to write the final ansatz parameters after Li simulation, in the param input format                                          |
| `out_truewavef`  | filename to which to cache the final true wavefunction (*not* that produced by the ansatz). This can be read by                                |
| `out_data`       | filename to which to write the simulation data, readable by Mathematica as an Association                                                      |
| [`in_truewavef`] | optional filename containing the true initial wavefunction. This is used to compare the simulations to true evolution when not starting at t=0 |

Simulation can be additionally configured by adjusting the below constants in `simulator.c`:

| constant               | explanation                                                                                                                                                                                                                                          |
|------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `TROT_CIRC_FN`         | filename of the input (single) Trotter cycle. This must correspond to the terms in the simulation Hamiltonian                                                                                                                                        |
| `SIM_HAMIL_FN`         | filename of the input simulation Hamiltonian                                                                                                                                                                                                         |
| `REAL_EVO_TIME_STEP`   | the timestep of the realtime simulation, i.e. how forward in time each iteration brings the simulated state. This controls the accuracy of the Li method, and the times the Trotter method will sample.                                              |
| `REAL_EVO_NUM_ITERS`   | the total number of iterations to perform. The total simulated duration is `REAL_EVO_TIME_STEP * REAL_EVO_NUM_ITERS`                                                                                                                                 |
| `REAL_EVO_TROT_MODE`   | controls how the Trotter method should be performed. `0` means every Trotter cycle is forward applied. `1` means every 2nd Trotter cycle is reversed. `2` means the total Trotter gate order is randomised. We find that `1` yields the best results |
| `REAL_EVO_TROT_CYCLES` | the number of Trotter cycles the Trotter method will use to produce a future state. `0` means the Trotter method isn't simulated (to speedup the code)                                                                                               |
| `RAND_TROT_SEED`       | seeds the randomness, only used when `REAL_EVO_TROT_MODE=2`                                                                                                                                                                                          |
| `OUTPUT_NUM_SIGFIGS`   | controls the number of digits in the mantissa of the output numbers written to the output files                                                                                                                                                      |

Furthermore, the following constants in `paramevolver.c` can be adjusted:

| constant                     | explanation                                                                                                                                                                                            |
|------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `REAL_SOLVER_METHOD`         | whether to use TSVD (`=0`) or Tikhonov regularisation (`=1`) to solve Li's variational equations each time-step. We recommend `1` to constrain that the parameters vary smoothly.                      |
| `TSVD_TOLERANCE`             | only applies if `REAL_SOLVER_METHOD=0`, in which case it controls what size singular values are neglected during TSVD. Singular values smaller than `TSVD_TOLERANCE * maxSingularValue` are discarded. |
| `TIKHONOV_PARAM_SEARCH_SIZE` | how many different Tikhonov regularisation parameters to test when choosing the optimal via the L-curve corner method                                                                                  |
| `DERIV_ORDER`                | what order (`1` to `4`) of finite difference formulae to use when approximating the wavefunction derivatives                                                                                           |
| `DERIV_STEP_SIZE`            | the parameter stepsize to use when approximating a wavefunction derivative by finite difference formulae                                                                                               |


### Output

After simulation, several files are created:

| filename        | explanation                                                                                                                             |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| `out_params`    | contains the final Li ansatz parameters, i.e. those which when fed into `ansatz`, approximate the wavefunction in file `out_truewavef`  |
| `out_truewavef` | contains amplitudes of the final true wavefunction `psi(t=REAL_EVO_TIME_STEP*REAL_EVO_NUM_ITERS)`                                       |
| `out_data`      | the main output file containing the simulation results, which can be read by Mathematica                                                |


`out_data` is actually an [Assocation](https://reference.wolfram.com/language/guide/Associations.html), which can be read in Mathematica...
```Mathematica
data = Get["out_data"]
````
and which contains the following keys (accessible by `data["key"]`, and reported by `Keys @ data`):

| key                      | explanation                                                                                                                                                                                  |
|--------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `variationalFidelityEvo` | fidelity of the wavefunctions produced by Li's method against the true evolution                                                                                                             |
| `trotterFidelityEvo`     | fidelity of the wavefunctions produced by Trotter's method against the true evolution                                                                                                        |
| `energyEvo`              | energy of the wavefunctions produced by Li's method according to `SIM_HAMIL_FN`. This is ideally constant when simulation is accurate                                                        |
| `paramsEvos`             | values of every parameter at each iteration of Li's method                                                                                                                                   |
| `trotterAngleEvos        | values of the parameters that Trotter's method assigns to the gates to produce each iteration's wavefunction                                                                                 |
| `residualEvo`            | the values of the residuals (norm for Tikhonov, sum of squares for TSVD) reported by GSL when solving the linear equations which update the parameters. This should grow as fidelity worsens |
| `hamilSpectrum`          | energy eigenvalues of `SIM_HAMIL_FN`                                                                                                                                                         |

------------------------------------------------------------------------------------------------


## Recompiler

### Usage

### How it works

### Output

------------------------------------------------------------------------------------------------

## Eliminator

### Usage

### How it works

### Output

------------------------------------------------------------------------------------------------

## Demonstration




