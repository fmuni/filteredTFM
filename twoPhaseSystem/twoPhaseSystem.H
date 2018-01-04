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

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Contributors
    Federico Municchi, TUGraz, 2017
\*---------------------------------------------------------------------------*/

#ifndef twoPhaseSystem_H
#define twoPhaseSystem_H

#include "IOdictionary.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dynamicParameters.H"
#include "auxEquations.H"
#include "stressClosure.H"
//#include "dragClosure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class heatTransferClosure;
class AuxEquations;
//class DynamicParameters;
//class StressClosure;
class dragClosure;


/*---------------------------------------------------------------------------*\
                      Class twoPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseSystem
:
    public IOdictionary
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to gravity
        const dimensionedVector&    g_;

        //- True to write internal fields
        bool   writeInternalFields_;

        //- Phase model 1
        phaseModel phase1_;

        //- Phase model 2
        phaseModel phase2_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- Dilatation term
        volScalarField dgdt_;

        //- Equation manager to solve additional equations
        autoPtr<AuxEquations> equationManager_;

        //- Advanced run-time calculations
        autoPtr<DynamicParameters> runTimeCalculations_;

        //- Optional dispersion diffusivity
        tmp<surfaceScalarField> pPrimeByA_;

        //- Phase pair for phase 1 dispersed in phase 2
        autoPtr<orderedPhasePair> pair1In2_;

        //- Drag closure
        autoPtr<dragClosure> drag_;

        //- Heat transfer model
        autoPtr<heatTransferClosure> heatTransfer_;

    // Private member functions

        //- Return the mixture flux
        tmp<surfaceScalarField> calcPhi() const;


public:

    // Constructors

        //- Construct from fvMesh
        twoPhaseSystem(const fvMesh&, const dimensionedVector& g);


    //- Destructor
    virtual ~twoPhaseSystem();


    // Member Functions

        //- Write relevant fields to disk
        void writeInternalFields()  const;

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Return the drag coefficient
        tmp<volTensorField> Kd() const;

        //- Return the face drag coefficient
        tmp<surfaceTensorField> Kdf() const;

        //- Return the average drag coefficient
        //- Each component is averaged over the corresponding
        //- UU tensor compenent
        tmp<volScalarField> KdAve() const;

        //- Apply drag force to the system of equations
        void applyDrag(volVectorField& U, fvVectorMatrix& UEqn);

        //- Return the face average drag coefficient
        //- Each component is averaged over the corresponding
        //- UU tensor compenent
        tmp<surfaceScalarField> KdfAve() const;

        //- Return the heat transfer coefficient
        tmp<volScalarField> Kh() const;

        //- Return sub-grid dispersion diffusivity
        tmp<volScalarField> D() const;

        //- Solve for the two-phase-fractions
        void solve();

        //- Correct two-phase properties other sub-grid stresses
        void correct();

        //- Correct two-phase turbulence
        void correctStressClosures();

        //- Read base phaseProperties dictionary
        bool read();

        // Access
            //- Reuturn calculator
            const DynamicParameters& dynPar() const
            {
                return runTimeCalculations_();
            }

             //- Reuturn calculator
             const AuxEquations& equations() const
             {
                 return equationManager_();
             }


             //- Reuturn calculator
             DynamicParameters& dynPar()
             {
                 return runTimeCalculations_();
             }

              //- Reuturn calculator
             AuxEquations& equations()
              {
                  return equationManager_();
              }


            //- Return the drag model for the given phase
            const dragClosure& drag(const phaseModel& phase) const;

            //- Return the surface tension coefficient
            const dimensionedScalar& sigma() const;

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the gravity
            inline const dimensionedVector& g() const;

            //- Return phase model 1
            inline const phaseModel& phase1() const;

            //- Return non-const access to phase model 1
            inline phaseModel& phase1();

            //- Return phase model 2
            inline const phaseModel& phase2() const;

            //- Return non-const access to phase model 2
            inline phaseModel& phase2();

            //- Return the phase not given as an argument
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Return the mixture flux
            inline const surfaceScalarField& phi() const;

            //- Return non-const access to the the mixture flux
            inline surfaceScalarField& phi();

            //- Return the dilatation term
            inline const volScalarField& dgdt() const;

            //- Return non-const access to the dilatation parameter
            inline volScalarField& dgdt();

            //- Return non-const access to the dispersion diffusivity
            inline tmp<surfaceScalarField>& pPrimeByA();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
