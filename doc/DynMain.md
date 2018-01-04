Dynamic Parameters
==

Description
--
This module allows the definition of specific parameters that are used used runtime.

Syntax
--

Dynamic parameters are defined in the _dynamicParameters_ sub-dictionary in the
[_phaseProperties_ dictionary](phasePropertiesDict.md).

* __constantFilterSize__ (optional): use this keyword when the mesh is uniform. When the keyword is present faster algorithms are used when filter size is required.

* __referenceFilterSize__: requires a number. Specifies the reference filter size to be used in some correlations.

* __settlingDrag__: requires a string. It's the drag law used to evaluate the settling velocity which is used in some correlations. Valid drag laws are: __Stokes__, __WenYu__, __KochHill__ and __Beetstra__.

Examples
--

```
dynamicParameters
{
  //Speed up calculations with constant cell size
  constantFilterSize;

  //Quantities used in correlations
  referenceFilterSize            1e-4;

  //Settling drag law for terminal velocity
  settlingDrag Stokes;
}

```
Back to [main index](01_main.md).
