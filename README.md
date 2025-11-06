Code for SuperDARN Borealis 3.5s integration into [lompe](https://github.com/klaundal/lompe), with plotting functions, otherwise known as the Fast Borealis Ionosphere (FBI). Can also be used to read in SuperDARN FitACF level data and put it into a format that Lompe accepts.

FBI data can be obtained from [superdarn.ca](https://superdarn.ca/fbi), which can be read in with this package.

# Installation
Install with the following command:

```
pip install git+https://github.com/billetd/FBI.git
```

All pre-requisites will be installed, if they are not already. 

**Important note:** Special forked branches of Lompe and pydarnio are also installed, with some small changes to make FBI work properly. If you already have `pydarnio` or `Lompe` installed, uninstall them before installing this.

# Creating an FBI output file

See [process_example.py](FBI/process_example.py) for an example of creating an FBI HDF5 file from SuperDARN data. The code was written to work with Borealis wide-beam data, but should work with any SuperDARN operating mode. 

The HDF5 file created contains the Lompe fits to the SuperDARN data used as input. See [readwrite.py](FBI/readwrite.py) for how this is formatted.

Code for making plots is also included. See [plotting_example.py](FBI/plotting_example.py) for an example of reading in an FBI HDF5 file and making plots.

# Basic usage for reading FitACF into a Lompe format

If all you want to do is read in fitacf SuperDARN files and incorporate them into lompe, along with data from other sources, your code should look something like this:
```python
import datetime as dt
import apexpy
import lompe
from FBI import fitacf
from FBI.process import prepare_lompe_inputs


all_data = fitacf.read_fitacfs(fitacfs, cores=5)  # This reads in data from a list of fitacf files you make
time = dt.datetime(2016, 5, 6, 2)  # Time of interest. Change accordingly.
apex = apexpy.Apex(time, refh=300)  # Make an apexpy object for coordinate transforms 
superdarn_data, _ = prepare_lompe_inputs(apex, all_data, time, 120, True)  # Make the lompe data object

# Pass to Lompe and make a fit. This assumes you’ve already defined conductances and have your own grid set up (See Lompe github for more details).
model = lompe.Emodel(grid, Hall_Pedersen_conductance=(SH, SP)) 
model.add_data(superdarn_data)  # Add other data objects as appropriate
model.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy') # Adjust regularisation as appropriate
```

