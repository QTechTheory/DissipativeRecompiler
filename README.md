DissipativeRecompiler
=====================

> **This README is a work in progress!**

**Table of Contents**
  - [Overview](#overview)
  - [Simulator](#simulator)
  - [Recompiler](#recompiler)
  - [Eliminator](#eliminator)
  - [Demonstration](#demonstration)

------------------------------------------------------------------------------------------------

## Overview

This repository contains C code for 3 executables.

- `simulator.c` which emulates real-time simulation using the [Li method](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.7.021050), Trotter method, and direct Hamiltonian exponentiation.
- `recompiler.c` which uses [imaginary-time simulation](https://arxiv.org/abs/1804.03023) to recompile one quantum circuit into another.
- `eliminator.c` which uses constrained imaginary time simulation to compactify a quantum circuit.

These can be compiled using [GNU Make](https://www.gnu.org/software/make/) with the 3 included makefiles, though this requires installing the GSL as outlined [here](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903).

The repository also includes a copy of [QuEST](https://quest.qtechtheory.org/) for simulating quantum circuits, in addition to

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


`out_data` is actually an [Assocation](https://reference.wolfram.com/language/guide/Associations.html), which can be read in Mathematica via
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
| `trotterAngleEvos`        | values of the parameters that Trotter's method assigns to the gates to produce each iteration's wavefunction                                                                                 |
| `residualEvo`            | the values of the residuals (norm for Tikhonov, sum of squares for TSVD) reported by GSL when solving the linear equations which update the parameters. This should grow as fidelity worsens |
| `hamilSpectrum`          | energy eigenvalues of `SIM_HAMIL_FN`                                                                                                                                                         |

------------------------------------------------------------------------------------------------


## Recompiler

`recompiler.c` takes a target wavefunction (produced by some parameters in some ansatz acting on some input state) and attempts to find parameters for a new ansatz which approximates that wavefunction. As input, we take a fictitious "recompilation" Hamiltonian which has the input state as its ground state, the old ansatz and its assigned parameters, and the new ansatz. The code then emulates the process of using variational imaginary time simulation to evolve the new ansatz parameters, such that the produced state approaches the target state. After many iterations, the new ansatz should reproduce the target state, and so be used in place of the input ansatz to generate that state. This is useful when the new ansatz is easier to implement (on quantum hardware) than the original ansatz, i.e. by being shorter, less sophisticated, or more robust to noise. It may even be useful for checkpointing dynamical simulation, by periodically compressing the flexible ansatz circuits into shorter ones only good for generating the current wavefunction, and extending them with more flexible gates (e.g. Trotter cycles).

`recompiler.c` can utilise an optional *luring* subroutine, whereby an intermediate states are targeted (by slowly activating the old ansatz parameters) before the ultimate target wavefunction. This can slow the convergence rate of the new parameters, but increases their robustness in reaching their optima.

### How it works

We accept an input state which is known to be the groundstate of our recompilation Hamiltonian. We also accept an old circuit which acts on the input state to produce the target state, which we desire a new ansatz to reproduce. Finding parameters of the new ansatz which, when acting on the input state, produce the target state is a problem we recast to one of energy minimisation: we seek new parameters which produce the lowest energy state in a modified circuit, namely the inverse of the new ansatz appended to the old ansatz.

We drive toward such a state using variational imaginary time evolution, which monotonically decreases the parameterised state's energy. Our numerical emulation of imaginary time evolution is nearly identical to that for realtime Li simulation. When energy has stabilised and is near to the known ground state, we can reverse the new ansatz and negate its final parameters to produce the target state from the input state. 

If convergence halted before the ground state was reached, we may employ luring. Before simulation, we scale down the old ansatz parameters, shifting that unitary "towards identity" and the state it produces toward the input state. This is an easier intermediate target state for the new ansatz to "reverse engineer" since it means the initial new state (identity acting on the old ansatz) is closer to ground. We can now perform imaginary time simulation and once reaching a sufficiently low energy, scale the old ansatz parameters toward their ultimate total values, raising the energy. We can repeat this to incrementally *lure* the new ansatz parameters toward their optima.

### Usage

After compiling with
```bash
make -f makefile_recompiler
```
the recompiler is run by suppling command line arguments
```bash
./recompiler old_ansatz old_params new_ansatz new_init_params true_state new_final_params out_data lures threshold
````
where

| argument      | explanation    |
|---------------|----------------|
|`old_ansatz`	| filename of the original ansatz circuit |
|`old_params`	| filename of the original ansatz parameters |
|`new_ansatz`	| filename of the new (recompiled) ansatz circuit |
|`new_init_params`	| filename of the inital values of the new ansatz circuit (typically ~0) |
|`true_state`	| filename of a wavefunction which will be compared to the state produced by the new ansatz, for fidelity logging purposes. This is typically the *true* future state of some realtime simulation, of which the old ansatz is an approximation |
| `new_final_params`	| filename to which to write the final new ansatz parameters. These are the parameters which, when utilised by the new ansatz, transform the input state into the target state (that produced by the old ansatz given the old parameters). This file is in the input parameter format |
| `out_data` | filename to which to write the simulation data, readable by Mathematica as an Association |
| `lures` | the number of intermediate *luring* states to target before targeting the ultimate old ansatz parameters. This can be `0` to perform no luring |
| `threshold` | only used when `lures > 0`, where it decides at what distance from ground state an intermediate luring stage is considered complete and the old parameters incremented toward their ultimate values |

There are some additional constants in `recompiler.c`:

| constant                     | explanation                                                                                                                                                                                                                             |
|------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `RECOMP_HAMIL_FN`            | filename of the input recompilation Hamiltonian whose groundstate is the input state to the old ansatz                                                                                                                                  |
| `IMAG_EVO_TIME_STEP`         | the timestep of the imaginary time simulation, performed with an imaginary adaptation of Li's algorithm                                                                                                                                 |
| `IMAG_EVO_NUM_ITERS`         | the number of iterations of variational imaginary time evolution to perform                                                                                                                                                             |
| `IMAG_EVO_INIT_TARGET_PARAM` | only used if `lures > 0` in which case it is the initial value of the old ansatz parameters. This must be kept very small but non-zero, since exactly zero initial parameters makes the first iteration's linear equations undetermined |

and a single constant in `paramevolver.c`: `IMAG_SOLVER_METHOD` which sets whether to use TSVD (`=0`) or Tikhonov regularisation (`=1`) when solving the linear equations in Li's algorithm each iteration. We find TSVD offers the most rapid and reliable convergence to the ground state.


### Output

After recompilation, several files are created:

| filename | explanation |
|----------|-------------|
| `new_final_params` | contains the final values of the new params after imaginary time simulation (in the input param format). These should be reversed and negated, and the new ansatz reversed, in order to apply them to reproduce the target state |
| `out_data` | the main output file containing the simulation results, which can be read by Mathematica |

`out_data` is actually an [Assocation](https://reference.wolfram.com/language/guide/Associations.html), which can be read in Mathematica via
```Mathematica
data = Get["out_data"]
````
and which contains the following keys (accessible by `data["key"]`, and reported by `Keys @ data`):

| key | explanation |
|-----|-------------|
| `newParamEvos` | values of every new parameter at each iteration of imaginary time simulation |
| `fidelityWithFinalEvo` | fidelity of the old ansatz state against the new ansatz state. This should increase with imaginary time iterations |
| `energyWithFinalEvo` | the energy of the state produced by applying the old ansatz then the inverse of the new ansatz. This should approach the groundstate of `RECMOP_HAMIL_FN` |
| `retargetIters` | populated if `lures>0`, in which case it is a list of iterations where the target state was adjusted because the current energy dropped below `threshold` |
| `fidelityWithTargetEvo` | populated if `lures>0`, where it's the fidelity between the new ansatz state and the current target state which gradually transforms from the input state to the old ansatz state |
| `energyWithTargetEvo` | populated if `lures>0`, where it's the energy of the state produced by applying the inverse of the new ansatz state to the target state. This is the energy which controls the lure |
| `fidelityWithTrueEvo` | the fidelity of the new ansatz state with `true_state` |
|`hamilSpectrum`| energy eigenvalues of `RECOMP_HAMIL_FN` |

------------------------------------------------------------------------------------------------

## Eliminator

`eliminator.c` can be performed after using `recompiler.c` to further shrink the new ansatz. It selectively changes some gate parameters to be zero at a small cost in the new ansatz's ability to reproduce the target state. Gates with zero parameters (which effect the identity) can then be removed from the circuit. The maximum acceptable loss of fidelity, as monitored by the energy increase, decides how many gates can be removed.

### How it works

`eliminator.c` works just like `recompiler.c`: variational imaginary time simulation drives the system toward the groundstate of a fictitious Hamiltonian, so that the new ansatz parameters approach values which reproduce the target state (that produced by the old ansatz). However, the parameter evolution equations have additional constraints.

We pick the smallest magnitude parameter in the new ansatz, and constrain the evolution equations such that this parameter approaches 0. This is as simple as reformulating *A x = b* to *A' x = b - col(A)* where a column of A (corresponding to the constrained parameter) has been populated with the chosen value and moved to the RHS. Our forced change in a parameter will increase the energy. Imaginary time simulation will then optimise the unconstrained parameters so as to return the system to ground state, or as close as can be reached. Instead of constraining the selected parameter change immediately to zero, we can bound its maximum change each iteration to slowly excite the attemptedly cooled state in an analog to *luring*.
Since we numerically constrain the selected parameter, we can force it to become precisely 0. Thereafter, we constrain that parameter to remain 0 for all future evolution, and choose the next smallest magnitude parameter to eliminate.
We can continue repeating this, eliminating gates one by one, until the energy rises above some pre-set threshold.

We choose to halt as soon as this energy threshold is reached, restore the parameters before that particular 'failed elimination', and exit. In theory, this energy check could be performed strictly after the parameter has become 0, or elimination re-attempted on another parameter, to eliminate additional gates. We expect this gives little reward since energy is seen to monotonically increase during gate elimination.

### Usage

After compiling with
```bash
make -f makefile_eliminator
```
the eliminator is run by suppling command line arguments
```bash
./eliminator old_ansatz old_params new_ansatz new_init_params true_state new_final_params out_data
````
where

| argument      | explanation    |
|---------------|----------------|
|`old_ansatz`	| filename of the original ansatz circuit |
|`old_params`	| filename of the original ansatz parameters |
|`new_ansatz`	| filename of the new (recompiled) ansatz circuit, to be further cmopressed |
|`new_init_params`	| filename of the current values of the new ansatz circuit which approximate the target state (old ansatz with old params). Some of these values are to become 0, and the rest likely modified |
|`true_state`	| filename of a wavefunction which will be compared to the state produced by the new ansatz, for fidelity logging purposes. This is typically the *true* future state of some realtime simulation, of which the old ansatz is an approximation |
| `new_final_params`	| filename to which to write the final new ansatz parameters, some of which will be exactly zero (having been eliminated). These are the parameters which, when utilised by the new ansatz, transform the input state into the target state (that produced by the old ansatz given the old parameters). This file is in the input parameter format |
| `out_data` | filename to which to write the simulation data, readable by Mathematica as an Association |

There are some additional constants in `eliminator.c`:

| constant                     | explanation                                                                                                                                                                                                                             |
|------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `RECOMP_HAMIL_FN`            | filename of the input recompilation Hamiltonian whose groundstate is the input state to the old ansatz                                                                                                                                  |
| `IMAG_EVO_TIME_STEP`         | the timestep of the imaginary time simulation, performed with an imaginary adaptation of Li's algorithm                                                                                                                                 |
| `IMAG_SOLVER_METHOD` | whether to use TSVD (`=0`) or Tikhonov regularisation (`=1`) when solving the linear equations in the imaginary-time adaptation Li's algorithm each iteration. We recommend TSVD. |
| `ELIM_PARAM_MAX` | the maximum magnitude a parameter may have and still be considered "eliminated". While this must be smaller than the initial values of the to-be-eliminated parameters, its specific value is of little importance since the parameters can be constrained to become exactly 0 (within machine epsilon) |
| `PARAM_CHANGE_MAX` | the maximum change a to-be-eliminated parameter may experience in a single iteration, to force the elimination process to happen slowly and smoothly. This means a parameter `p` takes `ceiling(abs(p) / PARAM_CHANGE_MAX)` iterations to be eliminated. |
| `ENERGY_INCREASE_TOL_FAC` | a multiple of the original gap between true ground state and the energy of the state produced by the new ansatz's initial parameters. When the energy reaches this multiplied by the original gap, the energy is considered to have risen unacceptably and the elimination subroutine is terminated. |
| `MAX_TOTAL_ITERS` | merely a bound on the total number of elimination iterations which can be performed, for memory purposes |

### Output

After elimination, two files are created:

| filename | explanation |
|----------|-------------|
| `new_final_params` | contains the final values of the new params after elimination, some of which are exactly 0. These should be reversed and negated, the new ansatz reversed and its 0-param gates removed, in order to reproduce the target state |
| `out_data` | the main output file containing the simulation results, which can be read by Mathematica |

`out_data` is actually an [Assocation](https://reference.wolfram.com/language/guide/Associations.html), which can be read in Mathematica via
```Mathematica
data = Get["out_data"]
````
and which contains the following keys (accessible by `data["key"]`, and reported by `Keys @ data`):

| key | explanation |
|-----|-------------|
| `origDist` | the distance of the original state (inverse new ansatz applied to old ansatz) from true ground |
| `totalIters` | the total number of imaginary-time iterations performed before the energy rose above `ENERGY_INCREASE_TOL_FAC * origDist - groundstate` and elimination was terminated |
| `totalParamChanges` | the total number of parameters eliminated, i.e. set exactly to 0 |
| `itersOfParamChanges` | the imaginary-time iterations after which a parameter was eliminated |
| `indsOfParamChanges` | the indices of the eliminated parameters (with respect to `new_init_params`) in order of their elimination time |
| `paramsEvo` | the value of every new parameter at every iteration |
| `energyEvo` | the energy of the parameterised state (inverse new ansatz applied to old ansatz) at every iteration |
| `fidelityInitEvo` | the fidelity between states created by the new ansatz and the old ansatz, at every iteration |
| `fidelityTrueEvo` | the fidelity between the state created by the new ansatz and `true_state` |

------------------------------------------------------------------------------------------------

## Demonstration




