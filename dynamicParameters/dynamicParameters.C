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
#include "dynamicParameters.H"
#include "twoPhaseSystem.H"
#include "phaseModel.H"
#include "error.H"

using namespace Foam;
//-------------------------- Constructors ---------------------------------//
Foam::DynamicParameters::DynamicParameters( const twoPhaseSystem&             fluid,
                                            dictionary   auxEquationsDict
                                          )
:
fluid_(fluid),
settings_(auxEquationsDict),
refDeltaf_(readScalar(settings_.lookup("referenceFilterSize"))),
phiPMax_(fluid_.phase1().alphaMax()),
filterWidth_(pow(fluid_.mesh().V()[0],1./3.)),
slipU_
{
     IOobject
    (
     "slipU",
     fluid_.mesh().time().timeName(),
     fluid_.mesh(),
     IOobject::NO_READ,
     IOobject::NO_WRITE
    ),
    fluid_.mesh(),
    dimensionedVector("0",dimless,vector::zero)
},
Sr_
{
     IOobject
    (
     "StrainRate",
     fluid_.mesh().time().timeName(),
     fluid_.mesh(),
     IOobject::NO_READ,
     IOobject::NO_WRITE
    ),
    fluid_.mesh(),
    dimensionedScalar("0",dimensionSet(0,0,-1,0,0,0,0),0.)
}
{
  //Set settling information
  word dragSettling(settings_.lookup("settlingDrag"));
  setupSettling(dragSettling, filterWidth_);

  if(settings_.found("constantFilterSize"))
  {
    FrFilter=&DynamicParameters::FrConstant;
    deltaF=&DynamicParameters::deltaFConstant;
  }
  else
  {
    FrFilter=&DynamicParameters::FrNoConstant;
    deltaF=&DynamicParameters::deltaFNoConstant;
  }

};
//-------------------------- Destructors ----------------------------------//
Foam::DynamicParameters::~DynamicParameters()
{
};
//---------------------------   Methods  ----------------------------------//
void Foam::DynamicParameters::setupSettling(word dragLaw, double filterSize )
{
  basicCalculations_.setupSettling(
                     fluid_.phase1().d()()[0],
                     fluid_.phase1().rho()()[0],
                     mag(fluid_.g()).value(),
                     fluid_.phase2().nu()()[0],
                     fluid_.phase2().rho()()[0],
                     dragLaw,
                     filterSize
  );
}

double Foam::DynamicParameters::FrConstant(int celli) const
{
  return basicCalculations_.settling.FrPf;
}

double Foam::DynamicParameters::FrNoConstant(int celli) const
{
  return (
            basicCalculations_.settling.FrPf
          * filterWidth_
          / pow(fluid_.mesh().V()[celli],1./3.)
         );
}

double Foam::DynamicParameters::filterFr(int celli) const
{
   return (this->*FrFilter)(celli);
}

double Foam::DynamicParameters::filterSize(int celli, int dim) const
{
   return (this->*deltaF)(celli,dim);
}

double Foam::DynamicParameters::deltaFConstant(int celli, int dim) const
{
 if(dim==DIMENSIONLESS)
  return filterWidth_/(settlingU()*settlingU())*gValue();

 return filterWidth_;
}

double Foam::DynamicParameters::deltaFNoConstant(int celli,int dim) const
{
 if(dim==DIMENSIONLESS)
  return pow(fluid_.mesh().V()[celli],1./3.)/(settlingU()*settlingU())*gValue();

 return pow(fluid_.mesh().V()[celli],1./3.);
}

const Foam::volScalarField& Foam::DynamicParameters::getProcessedField(
                                                         volScalarField& baseField,
                                                         word operationName
) const
{
  if(operationName == "none")
   return baseField;

  for(unsigned int id=0;id<filteredFieldsScalar_.size();id++)
  {
    if( (baseField.name()+ "." + operationName) == filteredFieldsScalar_[id]->name())
     return (*filteredFieldsScalar_[id]);
  }

  FatalError << "No operation " << operationName
             << " was defined for field " << baseField.name()
             << exit(FatalError);
  return baseField;
}

const Foam::volScalarField& Foam::DynamicParameters::getProcessedField(
                                                         word baseName,
                                                         word operationName
) const
{

  for(unsigned int id=0;id<filteredFieldsScalar_.size();id++)
  {
    if( (baseName + "." + operationName) == filteredFieldsScalar_[id]->name())
     return (*filteredFieldsScalar_[id]);
  }

  FatalError << "No operation " << operationName
             << " was defined for field " << baseName
             << exit(FatalError);
  //Return generic field
  return fluid_.phase1().rho();
}

void Foam::DynamicParameters::createOperationField(
                                                         volScalarField& baseField,
                                                         word operationName
) const
{
  if(operationName=="none")
   return;

  FatalError << "Could not find operation " << operationName << exit(FatalError);
}

double Foam::DynamicParameters::settlingU() const
{
  return basicCalculations_.settling.u;
}

const Foam::volVectorField&
Foam::DynamicParameters::scaledSlipU( const phaseModel& continuous,
                                      const phaseModel& dispersed
                                    ) const
{
    //TODO: the update should be in a separate function!
    dimensionedScalar dimSettling("set",
                                   dimensionSet(0,1,-1,0,0,0,0),
                                   settlingU()
                                 );

    slipU_ =   (continuous.U()  -  dispersed.U())
                    /  dimSettling;

    return slipU_;
}

const Foam::volScalarField&
Foam::DynamicParameters::strainRate( const phaseModel& phase
                                    ) const
{
    //TODO: the update should be in a separate function!
    Sr_ = sqrt(2.) * mag(symm(fvc::grad(phase.U())));
    return  Sr_;
}

double Foam::DynamicParameters::gValue() const
{
  return mag(fluid_.g()).value();
}

Foam::dimensionedVector Foam::DynamicParameters::g()  const
{
  return fluid_.g();
}

void  Foam::DynamicParameters::checkFilterSize(scalar refFilterSize) const
{
  double toleranceValue = 1e-10;

  forAll(fluid_.mesh().C(),celli)
  {
   double deltaF = double( filterSize(celli,DIMENSIONLESS) - refFilterSize );

    if(deltaF < toleranceValue)
    {

      FatalError << "\nA closure requires a reference filter size that is smaller"
               << " than the mesh size."
               << "\nYour mesh is " << filterSize(celli,DIMENSIONLESS)/refFilterSize
               << " times the reference filter size."
               << exit(FatalError);


   }
  }
}
