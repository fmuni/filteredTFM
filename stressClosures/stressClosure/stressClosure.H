/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Contributors
    Federico Municchi, TUGraz, 2017
\*---------------------------------------------------------------------------*/

#ifndef Stress_CLOSURE_H
#define Stress_CLOSURE_H

#include "fvCFD.H"
#include "phaseModel.H"
#include "stressSubClosure.H"
#include <vector>

namespace Foam
{

 class StressClosure
 {
  protected:

  phaseModel&                  phase_;

  dictionary                settings_;

  //- Meso scale stress closure
  autoPtr<StressSubClosure>       mesoScaleStressClosure_;

  //- Frictional stress closure
  autoPtr<StressSubClosure>      frictionalStressClosure_;

  //- Kinetic stress closure
  autoPtr<StressSubClosure>      microScaleStressClosure_;

  //- True if dispersed phase
  bool        isDispersed_;



  public:


    //Members

    //Methods

     //Construct from phaseModel object and dictionary
     StressClosure
     (
      phaseModel&              phase,
      dictionary         closureDict,
      bool               isDispersed
     );

     //Destructor
     ~StressClosure();

     //Correct sub closures
     void correct();

     //Write sub closure quantities
     void writeFields() const;


     //Return sub-grid kinetic energy
     tmp<volScalarField>   k() const;

     //Return sub-grid dynamic viscosity
     tmp<volScalarField>   muEff() const;

     //- Return the effective viscosity
     tmp<volScalarField> nuEff() const;

     //- Return the turbulence kinetic energy dissipation rate
     tmp<volScalarField> epsilon() const;

     //- Return the Reynolds stress tensor
     tmp<volSymmTensorField> R() const;

     //- Return the phase-pressure'
     // (derivative of phase-pressure w.r.t. phase-fraction)
     tmp<volScalarField> pPrime() const;

     //- Return the effective stress tensor
     tmp<volSymmTensorField> devRhoReff() const;

     //- Return the source term for the momentum equation
     tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

 };
}

#endif
