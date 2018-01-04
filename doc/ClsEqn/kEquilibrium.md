Kinetic energy equilibrium closure
===

Description
--

Provides algebraic equation for the equilibrium formulation of the kinetic energy
for the continuous or dispersed phases. Therefore, it closes a __scalar field__.

Syntax
--
* __phase__: requires a string. If __continuous__ it closes the continuous phase equation, otherwise if __dispersed__ it closes the dispersed phase equation.

* __otherPhaseK__: requires a string. This model requires both the continuous and dispersed phase kinetic energy equations to be solved. Therefore, the user should
enter the name of the kinetic energy field corresponding to the other phase.

__NOTE__: the user has to provide an equation for the __otherPhaseK__ field! (see example below).

Examples
--

In the [auxiliary equations sub-dictionary](EqnNew.md):
```

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

  equation scalarEquilibrium;

  subClosures  kEquilibrium;

  kEquilibrium
  {
    phase         continuous;
    otherPhasek   k.particles;
  }
}


```

Back to [equation closures](../ClsEqn.md).
