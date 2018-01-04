eulerianFilteredTFM
==
A filtered Two Fluid Model library and solver for OpenFOAM®

eulerianFilteredTFM is part of CFDEMcoupling, and has been developed within the 
[NanoSim Project](http://sintef.no/NanoSim) project by DCS Computing GmbH (Linz, Austria) and 
Graz University of Technology (F. Municchi, S. Radl).

Copyright Notice
------------------

- Copyright 2017 - Graz University of Technology (F. Municchi, S. Radl).

All rights reserved.

Note that some modules (or files) within the eulerianFilteredTFM package are not under the copyright of the above copyright holders. These modules (or files) are included in the eulerianFilteredTFM package with a valid open-source license. Please see the text in the source files, or the corresponding sub-directory, for more detailed information on copyright and the license that applies for these files.

License
-----------------
See the [LICENSE.md](LICENSE.md) file for details.

Warranty
-----------------
eulerianFilteredTFM is distributed in the hope that it will be applied and further developed by researchers and engineerings in academia or industry. However, eulerianFilteredTFM comes without any warranty, without even the implied warranty of merchantability or fitness for a particular purpose. 


Scope
---------------------------------------
__eulerianFilteredTFM__ is a C++ library based on the open source finite
 volume library OpenFOAM® developed to solve the three dimensional unsteady spatially averaged Euler-Euler equations that describe a two phase system composed of
 a continuous phase (fluid) and a dispersed phase (solid particles).

The library features several closures for the interphase transfer coefficients, the
phase stress tensor and additional user defined equations.

Important Hints for Usage
-----------------
- eulerianFilteredTFM is integrated in a CFDEMcoupling package. You should find appropriate solvers and test cases in the 'applications/solvers/' and 'tutorials' folders of your CFDEMcoupling package.
- compilation works similar as for CFDEMcoupling. For installation instructions, see 'doc/INSTALL.md'
- a detailed documentation of the whole eulerianFilteredTFM package is available in the 'doc' subfolder.
