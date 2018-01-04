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

#include "SchneiderbauerPirker.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "dynamicParameters.H"

#define TOL 1e-10
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HDragCorrections
{
    defineTypeNameAndDebug(SchneiderbauerPirkerCorrection, 0);
    addToRunTimeSelectionTable(HDragCorrection, SchneiderbauerPirkerCorrection, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HDragCorrections::SchneiderbauerPirkerCorrection::SchneiderbauerPirkerCorrection
(
    const dictionary& dict,
    const phasePair& pair
)
:
    HDragCorrection(dict, pair),
    settings_(dict),
    HSchneiderbauerPirker_
    (
      IOobject
      (
          "HSchneiderbauerPirker",
           alphaD().fluid().mesh().time().timeName(),
           alphaD().fluid().mesh()
      ),
      alphaD().fluid().mesh(),
      dimensionedTensor("zero", dimless, tensor::zero)
    ),
    operationAlpha(settings_.lookupOrDefault<word>("phiOperation","none"))

{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HDragCorrections::SchneiderbauerPirkerCorrection::~SchneiderbauerPirkerCorrection()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::HDragCorrections::SchneiderbauerPirkerCorrection::Hf() const
{

    const fvMesh&  mesh = alphaD().fluid().mesh();

    const volVectorField& uslip( markers().scaledSlipU(alphaC(),alphaD()) );


    forAll(mesh.C(),celli)
    {
        //Get filterSize
        scalar Df = markers().filterSize(celli);

        //Calculate f_phi
        scalar f_phi(alphaD()[celli]*(alphaD().alphaMax() +TOL -alphaD()[celli]));

        //Calculate f_delta
        scalar f_delta( pow(Df,5.)/(pow(Df,5.) + 1024.));

        scalar phiStar((alphaD()[celli]+TOL)/alphaD().alphaMax());
        //Calculate h_phi
        scalar h_phiden(  phiStar*phiStar*phiStar - 7.247*phiStar*phiStar
                         + 6.289*phiStar +0.384 );

        if(mag(h_phiden)<TOL)
             h_phiden=TOL;

        scalar h_phi(
                      (-6.743*phiStar*phiStar + 6.728*phiStar)
                     /h_phiden
        );

        //Calculate h_delta
        scalar h_deltaden(pow(Df,2.664)-25.89);

        if(mag( h_deltaden)<TOL)
         h_deltaden=TOL;

        scalar h_delta(pow(Df,2.664)/h_deltaden);

        //Update drag correction
        HSchneiderbauerPirker_[celli].xx() =
                                (1.- h_phi*h_delta)
                              * pow(mag(uslip[celli].x())+TOL,-9.*f_phi*f_delta);

        HSchneiderbauerPirker_[celli].yy() =
                                 (1.- h_phi*h_delta)
                               * pow(mag(uslip[celli].y()+TOL),-9.*f_phi*f_delta);

        HSchneiderbauerPirker_[celli].zz() =
                                  (1.- h_phi*h_delta)
                                * pow(mag(uslip[celli].z())+TOL,-9.*f_phi*f_delta);

    }

    return HSchneiderbauerPirker_;


}



// ************************************************************************* //
