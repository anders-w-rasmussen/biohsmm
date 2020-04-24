# peakmaker
Python package for making custom genome annotation algorithms with Hidden semi-Markov Models

# Installing

```bash
> pip install git https://github.com/anders-w-rasmussen/peakmaker
```

Now install eigen3 (eigen.tuxfamily.org), version 3+ (gonna fix this in the future)
```bash
> export EIGENPATH=/path/to/eigen3/ 
```
Next go to the folder where peak-maker is installed cd into directory util/
```bash
> g++ libbwdbwd.cpp -o libfwdbwd.so --shared -fPIC -DNDEBUG -O3 -I$EIGENPATH 
```

# Using the package

```python
import numpy as np
from peakmaker.model import 

print("hello world")
```




