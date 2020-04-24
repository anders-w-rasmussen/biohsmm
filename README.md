
# peakmaker
<img src=/doc_files/pyclassifier.png alt="drawing" width="275"/>

Python package for making custom genome annotation algorithms with Hidden semi-Markov Models



# Installing

```bash
> pip install git https://github.com/anders-w-rasmussen/peakmaker
```

You need to have eigen3 installed. (eigen.tuxfamily.org), version 3+ (gonna fix this in the future). 
Next go to the folder where peak-maker is installed cd into directory util/

```bash
> export EIGENPATH=/path/to/eigen3/ 
> g++ libbwdbwd.cpp -o libfwdbwd.so --shared -fPIC -DNDEBUG -O3 -I$EIGENPATH 
```

# Using the package

Let's start with simulating some data (pretend these are genomic signals of somekind)

```python
import numpy as np
from matplotlib import pyplot as plt

# Let's simulate some observations
region_length = 5000

# Fix random state
np.random.seed(19680801)

# Allow normal obs to be observations that are normally distributed perhaps
# some measurement of regulatory function (see CRISPR-decryptr)

normal_observations = np.random.normal(0, 1, region_length)
normal_observations[1000:1100] = np.random.normal(0.5, 1, 100)
normal_observations[4400:4600] = np.random.normal(-0.5, 1, 200)

# Allow binary observations to be observations of some active
# enhancer / silencer mark

binary_observations = np.zeros(region_length)
binary_observations[np.random.randint(1000, 1100, 10)] = 1
binary_observations[np.random.randint(4400, 4600, 20)] = 1

# Add some background frequency of seeing our histone mark
binary_observations[np.random.randint(0, region_length, 10)] = 1

# Assume some data is missing (replae with NANs)
normal_observations[np.random.randint(0, region_length, 2500)] = np.nan
binary_observations[np.random.randint(0, region_length, 2500)] = np.nan

# Plot these observations
x = np.arange(0, region_length)
plt.subplot(211)
plt.title('Normal Observations')
plt.plot(np.argwhere(np.isfinite(normal_observations == True)), normal_observations[np.argwhere(np.isfinite(normal_observations == True))])
plt.subplot(212)
plt.title('Bernoulli Observations')
plt.plot(np.argwhere(np.isfinite(binary_observations == True)), binary_observations[np.argwhere(np.isfinite(binary_observations == True))], c='r')
plt.show()

```




