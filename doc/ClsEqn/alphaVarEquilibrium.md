Phase volume fraction variance equilibrium closure
===

Description
--

Provides algebraic equation for the equilibrium formulation of the phase volume fraction variance
for the dispersed phases. Therefore, it closes a __scalar field__.

__NOTE__: despite the closure is applied to the dispersed phase, the variance
for the continuous phase is identical.

Syntax
--

* __dispersedPhaseK__: requires a string. That is the name of the granular temperature field for the dispersed phase.

__NOTE__: the user has to provide an equation for the __dispersedPhaseK__ field! (see example below).

Examples
--

In the [auxiliary equations sub-dictionary](EqnNew.md):
```

fieldsToSolve (k.particles, k.air, alphaVar);

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

  equation scalarEquilibrium;

  subClosures  kEquilibrium;

  kEquilibrium
  {
    phase         continuous;
    otherPhasek   k.particles;
  }
}

alphaVar
{

  equation scalarEquilibrium;

  subClosures  alphaVarEquilibrium;

  alphaVarEquilibrium
  {
    dispersedPhaseK k.particles;
  }
}


```

Back to [equation closures](../ClsEqn.md).
