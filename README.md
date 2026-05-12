# Landslides-Bennu

Finite-volume solver for shallow granular flow on a rotating, top-shaped asteroid
(modelled as a double cone). Simulates successive landslide (avalanche) events and
tracks the resulting change in spin rate through angular momentum conservation.

For the underlying theory see the reference paper included in the repository:
`KumarJFM_granular_flow_rotating_gravitating.pdf`

---

## Repository structure

| File | Description |
|---|---|
| `config.h` | Compile-time grid resolution (`RES = 200`) and `PI` |
| `bennu_input.txt` | All runtime parameters (editable without recompilation) |
| `bennu_solver.cpp` | Main solver — NOC central scheme, TVD-RK2 time integration |
| `gravity_field.ipynb` | Computes gravity field via stacked-disc integration; writes `grav45.txt` and `eta45.txt` |
| `post_processing.ipynb` | Reads solver output and produces plots of h, u, v, omega, mass |
| `impact_energy.ipynb` | Impact energy estimates and core erosion calculation |
| `grav45.txt` | Pre-generated normalised gravity magnitudes (200 values, 45° apex) |
| `eta45.txt` | Pre-generated gravity vector angles in radians (200 values, 45° apex) |

---

## Dependencies

### C++ solver

| Tool | Tested version |
|---|---|
| `g++` | 13.3.0 |
| C++ standard | C++17 (`-std=c++17`) |

### Python notebooks

| Package | Tested version |
|---|---|
| Python | 3.12.3 |
| numpy | 1.26.4 |
| scipy | 1.11.4 |
| mpmath | 1.2.1 |
| matplotlib | 3.6.3 |
| jupyter / nbconvert | any recent |

Install with:
```
pip install numpy scipy mpmath matplotlib notebook
```

---

## Reproducing results

Run the steps below **in order**. Steps 1 and 2 only need to be repeated if you
change the cone geometry (`apex_half_angle_deg`) or need higher-accuracy gravity.

### Step 1 — Generate the gravity field

Open and run all cells of `gravity_field.ipynb`.

This integrates the gravitational potential of the double-cone body over stacked
discs using elliptic integrals (ellipK, ellipE, ellipPi). It writes:
- `grav45.txt` — normalised gravity magnitude at each surface grid point
- `eta45.txt` — angle of the gravity vector from the axial direction (radians)

> **Note:** Pre-computed files are included in the repository. You can skip this
> step unless you change `apex_half_angle_deg` or `n_surface` in the notebook.
> `n_surface` in the notebook must equal `RES` in `config.h`.
>
> Runtime with default settings (`n_surface = 200`, `n_integ = 100`): ~1 minute.

### Step 2 — Compile the solver

```bash
g++ -std=c++17 -O2 -o bennu_solver bennu_solver.cpp
```

> To change the grid resolution, edit `RES` in `config.h` and recompile.
> All other parameters (friction angle, spin rate, number of avalanche events, etc.)
> can be changed in `bennu_input.txt` without recompiling.

### Step 3 — Run the solver

```bash
./bennu_solver
```

Reads `bennu_input.txt` and the gravity field files, then runs 25 successive
avalanche events. Terminal output:

```
t=0.0012375000  step=1  alpha=0.0000001240  omega=1.5000000000  dt=0.0012375000
...
t=250.0059375003  step=202025  alpha=-0.0000831444  omega=1.1445657849  dt=0.0012375000
```

Output files written (prefix set by `output_prefix` in `bennu_input.txt`):

| File | Contents |
|---|---|
| `run1_h.txt` | Flow height `h` — all cells, all output steps |
| `run1_u.txt` | Downslope velocity `u` |
| `run1_v.txt` | Azimuthal velocity `v` |
| `run1_omega.txt` | Rotation rate at each output step |
| `run1_time.txt` | Cumulative time at each output step |
| `run1_mass.txt` | Total regolith mass at each output step |
| `run1_ke.txt` | Kinetic energy at each output step |

> With default settings (`RES = 200`, `t_final = 10`, `n_avalanches = 25`,
> `output_interval = 50`): ~202 000 time steps, ~4 000 output snapshots,
> ~2 minutes on a modern laptop.

### Step 4 — Post-process and plot

Open and run all cells of `post_processing.ipynb`.

Set `run_prefix` and `RES` at the top of the notebook to match the values used
in `bennu_input.txt` and `config.h`. Produces inline plots of:
- Height profiles at selected times
- Final h, u, v profiles
- Spin rate and mass evolution
- Kinetic energy evolution

---

## Key parameters (`bennu_input.txt`)

| Parameter | Default | Description |
|---|---|---|
| `apex_half_angle_deg` | 45.0 | Cone half-apex angle ξ (degrees) |
| `friction_angle_deg` | 20.0 | Coulomb friction angle δ (degrees) |
| `omega_init` | 1.5 | Initial non-dimensional spin rate |
| `aspect_ratio` | 0.005 | Flow shallowness parameter ε |
| `regolith_depth` | 1.0 | Initial uniform layer thickness (non-dim) |
| `n_avalanches` | 25 | Number of successive landslide events |
| `t_final` | 10.0 | Simulation time per event (non-dim) |
| `output_interval` | 50 | Write output every N time steps |
| `output_prefix` | run1 | Prefix for all output files |

---

## Physical model summary

The solver integrates the depth-averaged shallow granular flow equations in
downslope coordinates on a cone surface. The conserved variables are

```
p = h·τ,   q = h·u·τ,   r = h·v·τ
```

where `h` is flow height, `u` downslope velocity, `v` azimuthal velocity, and
`τ = x·sin(ζ)` the radial coordinate. Source terms include slope gravity,
centrifugal and Coriolis accelerations, and Coulomb basal friction. Angular
momentum of the system (regolith + solid core) is conserved at each time step
to update the spin rate `ω`.

Numerical scheme: Kurganov-Tadmor non-oscillatory central flux with MUSCL
reconstruction (minmod limiter) and two-stage TVD Runge-Kutta time integration.
