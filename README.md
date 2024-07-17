Code for SuperDARN Borealis 3.5s integration into [lompe](https://github.com/klaundal/lompe), with plotting functions.

Install with the following command:

```
pip install git+https://github.com/billetd/FBI.git
```

All pre-requisites will be installed, if they are not already. Special forked branches of Lompe and pydarnio are also installed, with some small changes, so uninstall the regular versions (if you have them) before installing this.

If all you want to do is read in fitacf SuperDARN files and incorporate them into lompe, your code should look something like this:
```python
import datetime as dt
import apexpy
import lompe
from FBI import fitacf
from FBI.process import prepare_lompe_inputs


all_data = fitacf.read_fitacfs(fitacfs, cores=5)  # this reads in data from a list of fitacf files you make
time = dt.datetime(2016, 5, 6, 2)  # Time of interest. Change accordingly.
apex = apexpy.Apex(time, refh=300)  # Make an apexpy object for coordinate transforms 
superdarn_data, _ = prepare_lompe_inputs(apex, all_data, time, 120, True)  # Make the lompe data object

model = lompe.Emodel(grid, Hall_Pedersen_conductance=(SH, SP))  # Assuming you’ve already defined conductances and have your own grid set up
model.add_data(superdarn_data)  # Add other data objects as appropriate
model.run_inversion(l1=10, l2=0.1, lapack_driver='gelsy')
```

BUG WARNING

There is currently a bug in Apexpy ([issue here](https://github.com/aburrell/apexpy/issues/134)) that stops it working with Numpy version 2.0. Untill this is fixed, you will need to downgrade Numpy is order to use FBI. After installing FBI:
```python 
pip uninstall numpy
``` 
and then 
```python 
pip install numpy==1.26.4
``` 
