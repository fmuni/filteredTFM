eulerianFilteredTFM - Documentation
==

A filtered Two Fluid Model library for OpenFOAM®
--

__eulerianFilteredTFM__ is C++ library based on the open source finite
 volume library OpenFOAM® developed to solve the three dimensional unsteady spatially averaged Euler-Euler equations that describe a two phase system composed of
 a continuous phase (fluid) and a dispersed phase (solid particles).

The library features several closures for the interphase transfer coefficients, the
phase stress tensor and additional user defined equations.

Basic structure
--

The __eulerianFilteredTFM__ library stems from the __twoPhaseEulerFoam__ solver included
in OpenFOAM®, but its structure is more in line with the physics and the terminology used in fluid-solid suspensions.
The library can be subdivided in modules as follows:

* __twoPhaseSystem__: Similarly to _twoPhaseEulerFoam_, it is the main module that provides access to relevant quantities or operations and manages the two phases. The module is also solving the phase volume fraction equations.

* __interphase closures__: Provides closure models for the filtered interphase transfer coefficients (i.e., drag and heat transfer coefficients). Interphase closures are calculated as functions of:
 * __micro closures__: They are the microscopic closures corresponding to the functional form of the interphase transfer coefficients for the _unfiltered_ model.
 * __heterogeneity corrections__: They are corrections based on the fact that highly
 heterogeneous structures (like clusters) that affect the interphase exchange may exist at a sub-grid level.


* __stress closures__: Provides closure models for the phase filtered stress tensor.
Such tensor is expressed as the sum of three tensors representing the multiscale nature of the non-linear term in the transport equations:

 * __meso-scale stress__: At the meso-scale, only homogenous structures are considered.
These homogeneous structures are triggered by two kind of instabilities: one inertial
due to the relative motion between the gas and particle phases, and one diffusive that arise from the dumping of particle velocity fluctuations due to collisions or interstitial fluid. These structures are namely clusters and streamers of particles.
 * __frictional stress (particle phase only)__: It is the stress arising by frictional effects due to high particle packing. It is generally associated with band-like heterogeneous structures arising at a sub-grid scale.
 * __micro-scale stress__: It is the stress induced by kinetic phenomena like binary collisions (both for particle or fluid phases). This stress can be related to the formation of clusters at a sub-grid scale.


* __equation closures__: Additional closures for user-defined transport equations like granular temperature or phase fraction variance.

* __auxiliary equations__: This module allows to define additional equations to be solved alongside with mass and momentum conservation.

* __dynamic parameters__: Allows dynamic calculation of quantities to simplify the implementation of advanced closures.

Additionally, the library features the __eulerianFilteredTFM solver__ and verification cases.

Main index
--

|  Basics                     | Closures                            |    Auxiliary Equations               | dynamic Parameters  |
|:---                         |:---                                 |:---                                  |:---                 |
| [Installation](INSTALL.md)  | [Stress Closures](ClsStress.md)     | [Defining a new equation](EqnNew.md) | [Usage](DynMain.md) |
| [Dictionary](phasePropertiesDict.md)       | [Interphase Closures](ClsInter.md)  |            
|[Solvers](solvers.md)    | [Equation Closures](ClsEqn.md)      |                 
|  [Tutorials](tutorials.md)   | [Equation Closures](ClsEqn.md)      |     |
