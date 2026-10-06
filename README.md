Code for SuperDARN Borealis 3.5s integration into [lompe](https://github.com/klaundal/lompe), with plotting functions, otherwise known as the Fast Borealis Ionosphere (FBI). Can also be used to read in SuperDARN FitACF level data and put it into a format that Lompe accepts.

FBI data can be obtained from [superdarn.ca](https://superdarn.ca/fbi), which can be read in with this package.

# Installation
FBI requires Python 3.11+ to work correctly (Scipy dependence conflict with pyDARN)

Install with the following command:

```
pip install git+https://github.com/billetd/FBI.git
```

All pre-requisites will be installed, if they are not already. apexpy is built from source, so a Fortran compiler (e.g. gfortran) is needed.

FBI uses the main branch of [lompe](https://github.com/klaundal/lompe), and pyDARN's `ehn/gate2geo_vectorisation` branch until [SuperDARN/pydarn#450](https://github.com/SuperDARN/pydarn/pull/450) is merged. It no longer needs geodarn. If you previously installed FBI with the forked versions of lompe or geodarn, reinstall it (or use a fresh virtual enviroment).

# Creating an FBI output file

See [process_example.py](examples/process_example.py) for an example of creating an FBI HDF5 file from SuperDARN data. The code was written to work with Borealis wide-beam data, but should work with any SuperDARN operating mode. 

To process whole days, `FBI.extras.process_dates()` finds the fitacf files in a directory of the form `fitacfs_root/YYYY/MM/`, and writes one HDF5 file per two hours into `output_dir/YYYY/MM/`, skipping any that already exist.

The HDF5 file created contains the Lompe fits to the SuperDARN data used as input, with one group per scan. Each group has:

| Datasets | Contents |
| --- | --- |
| `v_e_model`, `v_n_model`, `mlats_model`, `mlons_model` | Fit velocity [m/s], magnetic east and north, on the lompe grid |
| `e_pot_model` | Electric potential [V] on the lompe grid |
| `v_e_darngrid`, `v_n_darngrid`, `mlats_darngrid`, `mlons_darngrid` | Fit velocity on the SuperDARN equal area grid |
| `v_e_los`, `v_n_los`, `mlats_los`, `mlons_los`, `rids` | The line-of-sight velocities that went into the fit, and the station id of each |
| `bound_mlats`, `bound_mlons` | Boundary of the fit |
| `scan_year` ... `scan_second`, `scan_millisec` | Time of the scan |

Velocities and potentials are compressed with HDF5's scale-offset filter, keeping them to within 0.5 m/s and 0.5 V. Read them back with `FBI.readwrite.fbi_load_hdf5()`. See [readwrite.py](FBI/readwrite.py) for more.

# Plotting

See [plotting_example.py](examples/plotting_example.py). `FBI.plotting.plot_main.plot_records()` plots a list of records from `fbi_load_hdf5()`, one image each, as `'vectors'`, `'potential'` or `'potential_polar'`. `FBI.extras.plot_fbi_files()` does a whole directory of FBI files, skipping what is already plotted, so it can be run daily.

# Parallel processing

Reading, processing and plotting take a `cores` argument for the number of worker processes. `None` uses every CPU available (or the `FBI_CORES` environment variable, if set), and `1` runs everything in the calling process, for debugging or profiling. Workers are forked, so this needs Linux or macOS, and scripts must run FBI from under an `if __name__ == '__main__':` guard.

Each worker uses one BLAS thread. Set `FBI_BLAS_THREADS` before importing FBI to change this, e.g. when running with `cores=1`.

# Basic usage for reading FitACF into a Lompe format

If all you want to do is read in fitacf SuperDARN files and incorporate them into lompe, along with data from other sources, your code should look something like this:
```python
import datetime as dt
import lompe
from FBI import fitacf
from FBI.process import prepare_lompe_inputs


all_data = fitacf.read_fitacfs(fitacfs, cores=5)  # This reads in data from a list of fitacf files you make
time = dt.datetime(2016, 5, 6, 2)  # Time of interest. Change accordingly.
superdarn_data, _ = prepare_lompe_inputs(all_data, time, 120, True)  # Make the lompe data object

# Pass to Lompe and make a fit. This assumes you’ve already defined conductances and have your own grid set up (See Lompe github for more details).
model = lompe.Emodel(grid, Hall_Pedersen_conductance=(SH, SP)) 
model.add_data(superdarn_data)  # Add other data objects as appropriate
model.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy') # Adjust regularisation as appropriate
```

# Tests

```
pip install -e ".[dev]"
pytest
```
