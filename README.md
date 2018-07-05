# quick-moebius

Simple repository for quick implementation of rendering and exploration of limit sets, tiling and orbits of groups and semigroups generated from a finite number of Moebius Transformations.

**Dependencies:**

* numpy
* matplotlib

**Example Usage:**

```python
import numpy as np
from mobius_maps import MoebiusTransform, SemiGroup

m_1 = MoebiusTransform((1+0j,0+1j,0+0j,-1+1j))
m_2 = MoebiusTransform(np.array([[0+1j,0+0j],[0+0j,0+1j]]))

sg = SemiGroup([m1,m2])
sg.compute_limits(3)
```
