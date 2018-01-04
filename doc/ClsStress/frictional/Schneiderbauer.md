Schneiderbauer frictional stress sub-closure
==

Description
--
This closure provides granular pressure and viscosity to model the frictional regime.
__For dispersed phase only__.

Syntax
--
Three parameters are required by this sub-closure:

* __I0__: requires a number.
* __muSt__: requires a number.
* __muC__: requires a number.

For a description of these parameters please consult the pdf documentation or the
appendix of Schneiderbauer et al. 2013 (DOI:10.1002/aic.14155).

Examples
--
In [stressClosure](../../ClsStress.md):

```
      frictional
      {
            type SchneiderbauerFrictional;

            I0      0.297;
            muSt    0.3819;
            muC     0.6435;
      }
```
Notice that we used the same values as Schneiderbauer et al. (2013) for the model parameters.

Back to [stress closures](../../ClsStress.md).
