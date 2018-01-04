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

#ifndef STRESS_SUBCLOSURE_H
#define STRESS_SUBCLOSURE_H

#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"
#include "dynamicParameters.H"

namespace Foam
{

 class StressSubClosure
 {
  protected:

   phaseModel&                phase_;

   dictionary              settings_;

   //-Closure viscosity
   autoPtr<volScalarField>        nuPrimePtr_;

   //-Closure phase-pressure
   autoPtr<volScalarField>         pPrimePtr_;

   //-Closure bulk viscosity
   autoPtr<volScalarField>    lambdaPrimePtr_;

   //-Anisotropic contribution to R (diagonal elements only)
   autoPtr<volVectorField>     aSigmaPtr_;

   //-Threshold below which the quantity is lineary smoothed to 0.
   scalar    phaseThreshold_;

   //-Allocate memory
   void createViscosity(word closureName);
   void createPressure(word closureName);
   void createLambda(word closureName);
   void createASigma(word closureName);

   void checkAutoPtr(word property,bool empty) const;

   //- Return the effective viscosity
   const volScalarField& nuPrime() const;

   //- Return the effective viscosity (non const)
   volScalarField& nuPrime();

   //- Return the bulk viscosity
   const volScalarField&  lambdaPrime() const;

   //- Return the bulk viscosity (non-const)
   volScalarField&  lambdaPrime();


   //- Return the bulk viscosity
   const volVectorField&  aSigma() const;

   //- Return the bulk viscosity (non-const)
   volVectorField&  aSigma();

   //-Return the phase-pressure
   const volScalarField& pPrime() const;

   //-Return the phase-pressure (non-const)
   volScalarField& pPrime();

   //Return surface phase phase-pressure
   tmp<surfaceScalarField> pPrimef() const;

   //Shorthand for markers
   const DynamicParameters& markers() const;

   //Necessary for template
   const phaseModel&   dispersedPhase_;

   //-Blends field value to get zero in the very dilute case
   template<class T>
   void smoothToZero(T& field) const
   {
    //Always based on the particle volume fraction (i.e., phase1)
    forAll(dispersedPhase_.U(),celli)
    {
       if(dispersedPhase_[celli] < phaseThreshold_)
       {
           //Linear smoothing
           field[celli] = dispersedPhase_[celli]/phaseThreshold_ * field[celli];

       }
    }

   }

  public:

       //- Runtime type information
    TypeName("StressSubClosure");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            StressSubClosure,
            dictionary,
            (
                const dictionary& dict,
                phaseModel& phase
            ),
            (dict,phase)
        );


    // Constructors

        StressSubClosure
        (
            const dictionary& dict,
            phaseModel& phase
        );

    // Destructors

       virtual ~StressSubClosure();

    // Selector

       static autoPtr<StressSubClosure> New
       (
           const dictionary& dict,
           phaseModel& phase,
           word subClosureType
       );

    //Methods

    //Write fields
    void writeFields() const;

    //Update total quantities in stressClosure
    void updateNu(volScalarField& nu) const;
    void updateP(volScalarField& p) const;
    void updatePf(surfaceScalarField& pf) const;
    void updateLambda(volScalarField& lambda) const;
    void updateASigma(volVectorField& ASigma) const;


    //Correct sub closure
    virtual void correct() {};

    //Return sub grid kinetic energy
    virtual  tmp<volScalarField>   k() const;

    //- Return the turbulence kinetic energy dissipation rate
    virtual tmp<volScalarField> epsilon() const;

    //- Return the Reynolds stress tensor
    virtual tmp<volSymmTensorField> R() const;
    virtual tmp<volSymmTensorField> R(volVectorField& U) const;

    //- Return the effective stress tensor
    virtual tmp<volSymmTensorField> devRhoReff() const;
    virtual tmp<volSymmTensorField> devRhoReff(volVectorField& U) const;

    //- Return the source term for the momentum equation
    virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;
 };

}

#endif
