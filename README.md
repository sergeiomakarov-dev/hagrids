# HaGrids

**HaGrids** is a Hamiltonian model for simplified island divertor configurations. It generates island-aligned computational grids and provides tools for magnetic field chaos studies. In real space, the modeled configuration represents island chains formed on top of circular flux surfaces. The generated field-aligned grids can be used as input for EMC3-EIRENE simulations.

## Reference

The model and its applications are described in:

> S. O. Makarov, F. Reimold, N. Maaziz, R. Davies, E. Rodríguez, V. R. Winters, and F. Onorati,
> **"Hamiltonian model for simplified island divertor grids and magnetic field chaos studies"**,
> *Physics of Plasmas* **33**, 072514 (2026).
> DOI: [10.1063/5.0334588](https://doi.org/10.1063/5.0334588)

If you use this code in your work, please cite the paper:

```bibtex
@article{Makarov2026HaGrids,
  author  = {Makarov, S. O. and Reimold, F. and Maaziz, N. and Davies, R. and Rodr{\'i}guez, E. and Winters, V. R. and Onorati, F.},
  title   = {Hamiltonian model for simplified island divertor grids and magnetic field chaos studies},
  journal = {Physics of Plasmas},
  volume  = {33},
  number  = {7},
  pages   = {072514},
  year    = {2026},
  doi     = {10.1063/5.0334588}
}
```

## Requirements

Python 3 with:

- `numpy`
- `scipy`
- `matplotlib`
- `mpmath`

Some auxiliary analysis scripts are also provided in MATLAB (`*.m`).

## Usage

### Island-aligned grid generation

1. Create `input_params.dat` (not tracked by git) based on `input_params_example.dat`.
2. Modify `input_params.dat` according to the required parameters (geometry, iota profile, island `m`/`n` numbers, perturbation amplitude, grid resolution, etc.).
3. Run:

   ```bash
   python island_alinged_grid.py
   ```

Example input files are provided in [input_params_examples/](input_params_examples/). Generated grids and data are written to the [output/](output/) directory.

### One-zone grid generation

For a one-zone grid (set `flag_type_init = 4` in the input parameters), run:

```bash
python one_zone_grid.py
```

## Repository overview

| File(s) | Purpose |
| --- | --- |
| `input_params.py`, `input_params_example.dat` | Input parameter parsing and example configuration |
| `fun_time_independent_mh.py`, `fun_time_dependent_mh.py` | Core routines for the time-independent and time-dependent magnetic Hamiltonian |
| `init_and_run.py`, `run_time_indepndent_magnetic_hamiltonian.py`, `run_time_depndent_magnetic_hamiltonian.py` | Field-line tracing / Hamiltonian integration drivers |
| `island_alinged_grid.py`, `one_zone_grid.py` | Island-aligned and one-zone grid generation |
| `create_points_in_core.py`, `create_points_in_island.py`, `create_points_in_pfr.py`, `create_points_one_zone.py` | Base-point generation in the core, island, and private-flux regions |
| `magnetic_feild_calculation.py` | Magnetic field evaluation |
| `check_boundary_points_island_alinged_grids.py`, `transform_theta_to_vartheta_island_alinged_grids.py` | Grid checking and coordinate transformation utilities |
| `plot_*.py`, `read_and_plot_*.py`, `fun_plot.py`, `write_poicare_to_text.py` | Plotting and Poincaré-section utilities |
| `*.m` | MATLAB counterparts of selected routines |
| `Time_independent_hamiltonian_pendulum.tex` | Notes on the time-independent Hamiltonian (pendulum) formulation |

## License

This project is licensed under the terms of the [LICENSE](LICENSE) file.
