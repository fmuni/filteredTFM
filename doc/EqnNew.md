Defining a new equation
==

__eulerianFilteredTFM__ allows the definition and solution of additional fields.
This allows, for example, to include kinetic energy or phase volume fraction variance
in the computation.

This is performed in the __auxiliary equation__ module.

Description
--

Allows definition of auxiliary equations to be used in further calculations.

Syntax
--

In [phaseProperties](phasePropertiesDict.md), additional equations are defined
using the _equations_ sub-dictionary.

Here, additional fields are created from:

* __fieldsToSolve__: requires a list of strings. Contains the names of all the additional fields sto solve.

Each field requires its own sub-dictionary where the related equation is specified
in :

* __equation__: requires a string. The equation name is composed of two words: the
__rank__ and the __structure__.

The __rank__ can be:

* __scalar__: the related field will be a scalar field.
* __vector__: the related field will be a vector field.
* __tensor__: the related field will be a tensor field of rank two.

While the structure:

* __Transport__: the equation as the structure of a transport equation with time derivative and convective term.
* __Equilibrium__: the equation as the structure of an algebraic equation.

__NOTE__: transport equations require an additional entry:

* __carrierPhase__: requries a string. It can be:
 * __continuous__: the carrier phase is the continuous phase.
 * __dispersed__: the carrier phase is the dispersed phase.
 * __mixture__: the field is transported by the mixture velocity field.

Each auxiliary equation requires the user to specify a [closure](ClsEqn.md).

Examples
--

```
equations
{

 fieldsToSolve (k.particles, k.air);

 k.particles
 {

  equation scalarEquilibrium;  

  subClosures  kEquilibrium;

  kEquilibrium
  {
    phase         dispersed;
    otherPhasek   k.air;
  }
 }

 k.air
 {

  equation scalarTransport;

  carrierPhase continuous;

  subClosures  ();
}


}


```

Back to [main index](01_main.md).
