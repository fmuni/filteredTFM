Installation
==

The library requires __OpenFOAM-4.1__ or __OpenFOAM-4.x__ to compile correctly. Other versions of OpenFOAM may necessitaty slight modification of the library and the solver.

The library is installed automatically in a CFDEMCoupling® environment using
__cfdemCompCFDEM__. Just ensure that the library name and solver are in the 'library-list.txt' and 'solver-list.txt' file of your CFDEMcoupling package (these files reside in 'src/lagrangian/cfdemParticle/etc/')

Otherwise, it is sufficient to type __wmake libso__ in the __eulerianFilteredTFM__ folder
to compile the library and __wmake__ in the solver folder to compile the solver.
