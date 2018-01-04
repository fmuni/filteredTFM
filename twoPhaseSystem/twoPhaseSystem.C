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

#include "twoPhaseSystem.H"
#include "heatTransferClosure.H"
#include "fvMatrix.H"
#include "surfaceInterpolate.H"
#include "MULES.H"
#include "subCycle.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcCurl.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fixedValueFvsPatchFields.H"
#include "HashPtrTable.H"
#include "dragClosure.H"
#include "stressClosure.H"
#include "auxEquations.H"
#include "dynamicParameters.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::twoPhaseSystem
(
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    g_(g),

    writeInternalFields_(
                           (
                             found("writeInternalFields")
                           )
                           ?
                           true : false
                         
                       ),

    phase1_
    (
        *this,
        *this,
         //Allow the user to change phase name
        {

         (
           found("setPhaseNames")
         )
         ?
         ( word(subDict("setPhaseNames").lookup("dispersed")) )
         :
         "dispersed"
        },
       true
    ),

    phase2_
    (
        *this,
        *this,
         //Allow the user to change phase name
        {

         (
           found("setPhaseNames")
         )
         ?
         ( word(subDict("setPhaseNames").lookup("continuous")) )
         :
         "continuous"
       },
       false
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->calcPhi()
    ),

    dgdt_
    (
        IOobject
        (
            "dgdt",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dgdt", dimless/dimTime, 0)
    )


{
    phase2_.volScalarField::operator=(scalar(1) - phase1_);



    pair1In2_.set
    (
        new orderedPhasePair
        (
            phase1_,
            phase2_,
            g
        )
    );


    equationManager_.set
    (
        new AuxEquations
         (
            *this,
            subDict("equations")
         )
    );

    runTimeCalculations_.set
    (
        new DynamicParameters
         (
            *this,
            subDict("dynamicParameters")
         )
    );


    // Closures

    drag_.set
    (
        new dragClosure
        (
            subDict("drag"),
            pair1In2_()
        )
    );


    heatTransfer_.set
    (
        new heatTransferClosure
        (
            subDict("heatTransfer"),
            pair1In2_()
        )
    );




}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::~twoPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::rho() const
{
    return phase1_*phase1_.thermo().rho() + phase2_*phase2_.thermo().rho();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::U() const
{
    return phase1_*phase1_.U() + phase2_*phase2_.U();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::calcPhi() const
{
    return
        fvc::interpolate(phase1_)*phase1_.phi()
      + fvc::interpolate(phase2_)*phase2_.phi();
}


Foam::tmp<Foam::volTensorField> Foam::twoPhaseSystem::Kd() const
{
    return drag_->K();
}


Foam::tmp<Foam::surfaceTensorField> Foam::twoPhaseSystem::Kdf() const
{
    return drag_->Kf();
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kh() const
{
    return heatTransfer_->K();
}



void Foam::twoPhaseSystem::solve()
{
    const Time& runTime = mesh_.time();

    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    const surfaceScalarField& phi1 = phase1_.phi();
    const surfaceScalarField& phi2 = phase2_.phi();

    const dictionary& alphaControls = mesh_.solverDict
    (
        alpha1.name()
    );

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    alpha1.correctBoundaryConditions();


    surfaceScalarField phic("phic", phi_);
    surfaceScalarField phir("phir", phi1 - phi2);

    tmp<surfaceScalarField> alpha1alpha2f;

    if (pPrimeByA_.valid())
    {
        alpha1alpha2f =
            fvc::interpolate(max(alpha1, scalar(0)))
           *fvc::interpolate(max(alpha2, scalar(0)));

        surfaceScalarField phiP
        (
            pPrimeByA_()*fvc::snGrad(alpha1, "bounded")*mesh_.magSf()
        );

        phir += phiP;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", dgdt_.dimensions(), 0.0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi_)*min(alpha1, scalar(1))
        );

        forAll(dgdt_, celli)
        {
            if (dgdt_[celli] > 0.0)
            {
                Sp[celli] -= dgdt_[celli]/max(1.0 - alpha1[celli], 1e-4);
                Su[celli] += dgdt_[celli]/max(1.0 - alpha1[celli], 1e-4);
            }
            else if (dgdt_[celli] < 0.0)
            {
                Sp[celli] += dgdt_[celli]/max(alpha1[celli], 1e-4);
            }
        }

        surfaceScalarField alphaPhic1
        (
            fvc::flux
            (
                phic,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        surfaceScalarField::Boundary& alphaPhic1Bf =
            alphaPhic1.boundaryFieldRef();

        // Ensure that the flux at inflow BCs is preserved
        forAll(alphaPhic1Bf, patchi)
        {
            fvsPatchScalarField& alphaPhic1p = alphaPhic1Bf[patchi];

            if (!alphaPhic1p.coupled())
            {
                const scalarField& phi1p = phi1.boundaryField()[patchi];
                const scalarField& alpha1p = alpha1.boundaryField()[patchi];

                forAll(alphaPhic1p, facei)
                {
                    if (phi1p[facei] < 0)
                    {
                        alphaPhic1p[facei] = alpha1p[facei]*phi1p[facei];
                    }
                }
            }
        }

        if (nAlphaSubCycles > 1)
        {
            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhic10(alphaPhic1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi_,
                    alphaPhic10,
                    (alphaSubCycle.index()*Sp)(),
                    (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                    phase1_.alphaMax(),
                    0
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase1_.alphaPhi() = alphaPhic10;
                }
                else
                {
                    phase1_.alphaPhi() += alphaPhic10;
                }
            }

            phase1_.alphaPhi() /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhic1,
                Sp,
                Su,
                phase1_.alphaMax(),
                0
            );

            phase1_.alphaPhi() = alphaPhic1;
        }

        if (pPrimeByA_.valid())
        {
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian(alpha1alpha2f()*pPrimeByA_(), alpha1, "bounded")
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            phase1_.alphaPhi() += alpha1Eqn.flux();
        }

        phase1_.alphaRhoPhi() =
            fvc::interpolate(phase1_.rho())*phase1_.alphaPhi();

        phase2_.alphaPhi() = phi_ - phase1_.alphaPhi();
        alpha1 = max(0.0, alpha1); //bound!
        alpha2 = scalar(1) - alpha1;
        phase2_.alphaRhoPhi() =
            fvc::interpolate(phase2_.rho())*phase2_.alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }
}


void Foam::twoPhaseSystem::correct()
{
    phase1_.correct();
    phase2_.correct();
}


void Foam::twoPhaseSystem::correctStressClosures()
{
    phase1_.stressClosure().correct();
    phase2_.stressClosure().correct();
}


bool Foam::twoPhaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        readOK &= phase1_.read(*this);
        readOK &= phase2_.read(*this);

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


const Foam::dragClosure& Foam::twoPhaseSystem::drag(const phaseModel& phase) const
{
    return drag_();
}


const Foam::dimensionedScalar& Foam::twoPhaseSystem::sigma() const
{
    return pair1In2_->sigma();
}

Foam::tmp<volScalarField> Foam::twoPhaseSystem::D() const
{
  //NOT IMPLEMENTED: returns zero
  return phase1_.stressClosure().pPrime()*0.0;
}

Foam::tmp<volScalarField> Foam::twoPhaseSystem::KdAve() const
{

   volVectorField    Umix(phase1_*phase1_.U() + phase2_*phase2_.U());
   volScalarField    UUmag( mag(Umix)*mag(Umix) );
   volScalarField    KdNotNorm( Kd()&&(Umix * Umix ));

   dimensionedScalar eps1("eps1",UUmag.dimensions(),1e-32);
   dimensionedScalar eps2("eps2",KdNotNorm.dimensions(),1e-32);


   return(
           (KdNotNorm + eps2)/( UUmag +eps1 )
   );



}

Foam::tmp<surfaceScalarField> Foam::twoPhaseSystem::KdfAve() const
{
 volVectorField    Umix(phase1_*phase1_.U() + phase2_*phase2_.U());
 volScalarField    UUmag( mag(Umix)*mag(Umix) );
 surfaceScalarField    KdNotNorm( Kdf()&&fvc::interpolate(Umix * Umix ));

 dimensionedScalar eps1("eps1",UUmag.dimensions(),1e-32);
 dimensionedScalar eps2("eps2",KdNotNorm.dimensions(),1e-32);


 return(
         (KdNotNorm + eps2)/fvc::interpolate( UUmag +eps1 )
 );

}

void Foam::twoPhaseSystem::applyDrag(volVectorField&   U,
                                     fvVectorMatrix& UEqn)
{
  //Isotropic component of drag tensor
  tensorField                fullDrag( Kd() );
  volScalarField         isoDrag(1.0/3.0*tr(Kd()) );

  const scalarField& V = mesh_.V();
  scalarField& Udiag = UEqn.diag();
  vectorField& Usource = UEqn.source();

  //Apply drag to matrix
  forAll(Udiag, celli)
  {
   //Isotropic drag is implicit
    Udiag[celli] += V[celli]*isoDrag[celli];
    //Anisotropic components are explicit
    Usource[celli] -= V[celli]*( (fullDrag[celli] - tensor::I*isoDrag[celli] )
                                  & U[celli]
                               );
 }

}

void Foam::twoPhaseSystem::writeInternalFields() const
{
  if(writeInternalFields_ )
  {
    drag_->writeFields();
    phase1_.stressClosure().writeFields();
    phase2_.stressClosure().writeFields();
  }
}

// ************************************************************************* //
