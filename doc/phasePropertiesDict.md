phaseProperties dictionary
==
Description
--
__phaseProperties__ is the main dictionary for using the library and share some similarities with the analogous in _twoPhaseEulerFoam_.

__NOTE__: _eulerianFilteredTFM_ also requires material properties to be specified in
_thermophysicalProperties_.

Syntax
--

By default, _eulerianFilteredTFM_ is labeling the two phases as __dispersed__ and __continuous__, but the user can adjust the names using:

* __setPhaseNames__: is a sub-dictionary that allows the user to rename the default phases as shown below.

It is possible to write internal fields (viscosities, phase pressures, drag, etc. ) declaring:

* __writeInternalFields__

```
setPhaseNames
{
    dispersed  particles;
    continuous air;
}

writeInternalFields;

```

For each phase, an additional sub-dictionary is required where the following entries are specified:

* __diameterModel__: only __constant__ is currently available.
* __constantCoeffs__: this sub-dictionary requires the specification of the particle diameter __d__.
* __alphaMax__: maximum allowable phase fraction (also used in some correlations).
* __residualAlpha__: tolerance value for the phase volume fraction.

Other entries are described in the [stress closure section](ClsStress.md).

```
particles
{
    //Particles are glass beads
    diameterModel constant;
    constantCoeffs
    {
        d               58.9e-6;
    }

    alphaMax        0.62;
    residualAlpha   1e-6;


    stressClosure
    {
        mesoScale
        {
            type  SarkarMeso;
            phase dispersed;
        }

        frictional
        {
            type SchneiderbauerFrictional;

            //Parameters for SchneiderbauerFrictional
            I0      0.297;
            muSt    0.3819;
            muC     0.6435;
        }

        microScale
        {
            type SarkarMicro;
        }
    }

}

Remaining entries are described in the documentation for [interphase closures](ClsInter.md), [auxiliary equations](EqnNew.md) and [dynamic parameters](DynMain.md).

```

Examples
--

A full _phaseProperties_ dictionary could be:

```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  4.x                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      phaseProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Rename phases with case-appropriate names
setPhaseNames
{
    dispersed  particles;
    continuous air;
}

particles
{
    //Particles are glass beads
    diameterModel constant;
    constantCoeffs
    {
        d               58.9e-6;
    }

    alphaMax        0.62;
    residualAlpha   1e-6;

    stressClosure
    {
        mesoScale
        {
            type  SarkarMeso;
            phase dispersed;
        }

        frictional
        {
            type SchneiderbauerFrictional;

            //Parameters for SchneiderbauerFrictional
            I0      0.297;
            muSt    0.3819;
            muC     0.6435;
        }

        microScale
        {
            type SarkarMicro;
        }
    }

}

air
{
    diameterModel constant;
    constantCoeffs
    {
        d               1;
    }

    residualAlpha   0;

    stressClosure
    {
        mesoScale
        {
            type  SarkarMeso;
            phase continuous;
        }

        microScale
        {
            type molecular;
        }
    }

}


drag
{
       microscopicDragLaw
       {
        type            GidaspowErgunWenYu;
        residualRe      1e-3;
       }

       heterogeneousCorrection
       {
        type        SchneiderbauerPirker;
       }

}


heatTransfer
{
      microscopicNusselt
      {
        type            RanzMarshall;
        residualAlpha   1e-3;
      }

      heterogeneousCorrection
      {
       type        none;
      }

}

equations
{

}

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
