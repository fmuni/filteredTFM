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

#include "CloeteCorrection.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "dynamicParameters.H"

#define TOL 1e-10
#define PI 3.14159265359
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HDragCorrections
{
    defineTypeNameAndDebug(CloeteCorrection, 0);
    addToRunTimeSelectionTable(HDragCorrection, CloeteCorrection, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HDragCorrections::CloeteCorrection::CloeteCorrection
(
    const dictionary& dict,
    const phasePair& pair
)
:
    HDragCorrection(dict, pair),
    settings_(dict),
    HCloete_
     (
        IOobject
        (
            "HCloete",
             pair_.dispersed().fluid().mesh().time().timeName(),
             pair_.dispersed().fluid().mesh()
        ),
        pair_.dispersed().fluid().mesh(),
        dimensionedTensor("zero", dimless, tensor::zero)
     ),
    operationAlpha(settings_.lookupOrDefault<word>("phiOperation","none")),
    x1_(38.530),
    x2_(40.940),
    x3_(1.2880),
    x4_(0.8498),
    x5_(0.2815),
    x6_(566.30),
    x7_(472.70),
    x8_(0.7652),
    x9_(0.4772),
    Dfr_(0.1286),
    y1_(36.360),
    y2_(16.380),
    y3_(0.6339),
    y4_(atan(0.6417)),
    y5_(0.6811),
    y6_(1.9360),
    y7_(0.7122),
    y8_(0.9235),
    y9_(0.1870)
{


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HDragCorrections::CloeteCorrection::~CloeteCorrection()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::HDragCorrections::CloeteCorrection::Hf() const
{
    //Slip velocity
    const volVectorField& uslip(markers().scaledSlipU(alphaC(),alphaD()));
    scalar TwooverPIcube((2./pow(PI,3.)));
    const dimensionedVector g(markers().g());
    dimensionedScalar magg(mag(g));

   //Check if filter size is compatible with this model
   markers().checkFilterSize(Dfr_);



    forAll(uslip, celli)
    {
      //Slip velocity parallel to gravity
      vector slipPV(
                     (
                       (uslip[celli]&g)
                       *
                       (
                        g/(magg*magg)
                       )
                     ).value()
                  );



      //Module of the slip velocity parallel to gravity
      scalar slipP( mag(slipPV) );

      //Module of the slip velocity normal to gravity
      scalar slipN( mag(uslip[celli]-slipPV));


      //DeltaFilter

      scalar DeltaF( markers().filterSize(celli)
                     -Dfr_
                   );



      scalar alphaStar( alphaD().alphaMax()
                      - alphaD()[celli]);
      //Correction parallel to gravity
      scalar cP(
                 exp(
                      -(
                           TwooverPIcube
                         * atan( x1_*pow(DeltaF,x8_)*alphaStar )
                         * atan( x2_*pow(DeltaF,x9_)*alphaStar)
                         *
                         (
                             x4_*log10(slipP+TOL) + x5_
                           + x6_*log10(slipP+TOL)*log10(slipP+TOL)
                           * (scalar(1.0) * 2.0/PI * atan(x7_*DeltaF))
                         )

                       )
                    )
      );

      //Correction normal to gravity
      scalar cN(
                 exp(
                      -(
                          TwooverPIcube
                         * atan( y1_*pow(DeltaF,y8_)*alphaStar )
                         * atan( y2_*pow(DeltaF,y9_)*alphaStar)
                         * atan( y3_*DeltaF)
                         *
                         (
                             y4_*slipN
                           / atan(y3_*DeltaF)
                           * x5_ *log10(slipN+TOL)
                           + x6_*alphaD()[celli]
                           + x7_
                         )

                       )
                    )
      );

      //Add to tensor diagonal
      vector comp(  cP*slipPV/(slipP+TOL)
                   + cN*(uslip[celli]-slipPV)/(slipN+TOL)
                  ) ;


      HCloete_[celli].xx()= 1.-mag(comp.x());
      HCloete_[celli].yy()= 1.-mag(comp.y());
      HCloete_[celli].zz()= 1.-mag(comp.z());

    }


    return HCloete_;
}



// ************************************************************************* //
