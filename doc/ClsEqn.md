Equation Closures
==

Description
--

Closures for user-defined equations introduced via the [auxiliary equations](EqnNew.md) utility. They are specified for each auxiliary equation.

Syntax
--

In _eulerianFilteredTFM_, each equation has one closure that is composed of several sub-closures. Such sub-closures are declared in the [equation sub-dictionary](EqnNew.md) using:

*  __subClosures__: requires a list of strings. List of sub closures used to close the auxiliary equation.

In addition, sub-closures are defined in dedicated sub-dictionaries.

Examples
--
In the [equation sub-dictionary](EqnNew.md):
```
 subClosures  kEquilibrium;

 kEquilibrium
 {
   phase         dispersed;
   otherPhasek   k.air;
 }
```
This example provides the equilibrium closure to the equilibrium closure for the granular temperature.

Available sub-closures
--

| scalar field | vector field | tensor field |
|:-- |:-- |:-- |
| [kinetic energy equilibrium](ClsEqn/kEquilibrium.md)
| [phase volume fraction variance equilibrium](ClsEqn/alphaVarEquilibrium.md)
