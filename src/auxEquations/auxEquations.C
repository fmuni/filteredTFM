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

#include "auxEquations.H"
#include "error.H"
#include "dynamicParameters.H"
#include "eqnClosure.H"


#define TOL 1e-20

using namespace Foam;
//-------------------------- Constructors ---------------------------------//
AuxEquations::AuxEquations( twoPhaseSystem&             fluid,
                            dictionary       auxEquationsDict
                           )
:
fluid_(fluid),
settings_(auxEquationsDict)
{
  //Fill information container
  if(settings_.found("fieldsToSolve"))
   eqnInfo_.names = wordList(settings_.lookup("fieldsToSolve"));

  forAll(eqnInfo_.names,eq)
  {
    dictionary eqnDict(settings_.subDict(eqnInfo_.names[eq]));

    word type(eqnDict.lookup("equation"));

    std::size_t found(0);

    //Check field rank and call sield constructor
    if(type.find("scalar")!=std::string::npos)
    {
       eqnInfo_.rank.push_back(SCALAR);
       found = type.find("scalar");
       volScalarField * tmp
       (
         new volScalarField
         (
              IOobject
              (
                  eqnInfo_.names[eq],
                  fluid_.mesh().time().timeName(),
                  fluid_.mesh(),
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              fluid_.mesh(),
              dimensionedScalar("0", dimless, 0. )
         )
       );

      scalarFields_.push_back(tmp);
    }
    else if(type.find("vector")!=std::string::npos)
    {
      eqnInfo_.rank.push_back(VECTOR);
      found = type.find("vector");

      volVectorField * tmp
      (
        new volVectorField
        (
             IOobject
             (
                 eqnInfo_.names[eq],
                 fluid_.mesh().time().timeName(),
                 fluid_.mesh(),
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             fluid_.mesh(),
             dimensionedVector("0", dimless, vector::zero )
        )
      );

     vectorFields_.push_back(tmp);
    }
    else if(type.find("tensor")!=std::string::npos)
    {
      eqnInfo_.rank.push_back(TENSOR);
      found = type.find("tensor");

      volTensorField * tmp
      (
        new volTensorField
        (
             IOobject
             (
                 eqnInfo_.names[eq],
                 fluid_.mesh().time().timeName(),
                 fluid_.mesh(),
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             fluid_.mesh(),
             dimensionedTensor("0", dimless, tensor::zero )
        )
      );

     tensorFields_.push_back(tmp);
    }
    else
    {
      FatalError << "Invalid equation " << type
                 << "for field" << eqnInfo_.names[eq];
    }

    //Check equation type
    if(type.find("Transport",found+1)!=std::string::npos)
    {
       eqnInfo_.type.push_back(TRANSPORT);

       //Assign transporting phase
       word trnsPhase(eqnDict.lookup("carrierPhase"));
       if(trnsPhase=="mixture")
       {
        eqnInfo_.phase.push_back(MIXTURE);
       }
       else if(trnsPhase=="continuous")
       {
        eqnInfo_.phase.push_back(CONTINUOUS);
       }
       else if(trnsPhase=="dispersed")
       {
        eqnInfo_.phase.push_back(DISPERSED);
       }
       else
       {
        FatalError << "Invalid equation " << type
                   << "for field" << eqnInfo_.names[eq] << "\n"
                   << "Carrier phase is not specified!";
       }

    }
    else if(type.find("Equilibrium",found+1)!=std::string::npos)
    {
      eqnInfo_.type.push_back(EQUILIBRIUM);

      eqnInfo_.phase.push_back(MIXTURE);
    }
    else
    {
      FatalError << "Invalid equation " << type
                 << "for field" << eqnInfo_.names[eq];
    }

    Info << "Found " << type << " equation for field "
         << eqnInfo_.names[eq];

    //Read closure
    EqnClosure * tmp = new EqnClosure(
                                       fluid_,
                                       settings_,
                                       eqnInfo_.names[eq]
                                     );

    eqnInfo_.closures.push_back(tmp);


  }


};
//-------------------------- Destructors ----------------------------------//
AuxEquations::~AuxEquations()
{
   for(unsigned int id=0;id<scalarFields_.size();id++)
    delete scalarFields_[id];

   for(unsigned int id=0;id<vectorFields_.size();id++)
    delete vectorFields_[id];

   for(unsigned int id=0;id<tensorFields_.size();id++)
    delete tensorFields_[id];

   for(unsigned int id=0;id<eqnInfo_.closures.size();id++)
    delete eqnInfo_.closures[id];


};
//---------------------------   Methods  ----------------------------------//
void AuxEquations::requestEquation(word equationName, word requestingClass) const
{
    forAll(eqnInfo_.names,eq)
    {
        if(equationName == eqnInfo_.names[eq])
         return;
    }

    FatalError << "Cannot find equation " << equationName
               << "requested by " << requestingClass << " !! ";
}

template<class T>
tmp<fvMatrix<T>> AuxEquations::transportEquation(
                                                  GeometricField<T,fvPatchField,volMesh>& field,
                                                  int eqId
                                                )

{

  if(eqnInfo_.phase[eqId] == MIXTURE)
  {
     //Just transported by the mixture velocity field
     return
     (
        fvm::ddt(field)
      + fvm::div(fluid_.phi(), field)
     );
  }

  phaseModel * alpha(NULL);

  if(eqnInfo_.phase[eqId] == CONTINUOUS)
  {
      //transported by the continuous phase
      alpha = &(fluid_.phase2());

  }
  else if(eqnInfo_.phase[eqId] == DISPERSED)
  {
      //transported by the dispersed phase
      alpha = &(fluid_.phase1());

  }

   //Return phase transport equation
    return
    (
       //Total derivative - main transport terms
       fvm::ddt((*alpha),alpha->rho(),field)
     + fvm::div(alpha->alphaRhoPhi(), field)
       //Additional term due to the chain rule
       //(i.e. derivatives of the alpha field)
       //are added explicitly
     - fvc::Sp(
                 fvc::ddt((*alpha), alpha->rho())
               + fvc::div(alpha->alphaRhoPhi())
               , field
               )

    );


}

template<class T>
tmp<fvMatrix<T>> AuxEquations::equilibriumEquation(
                                                    GeometricField<T,fvPatchField,volMesh>& field,
                                                    int eqId
                                                  )
{
    //Just returns an empty matrix
     return
    (
      fvMatrix<T>(field,dimless)
    );
}

void AuxEquations::solveField(word fieldName)
{

   forAll(eqnInfo_.names,eq)
   {
       if(eqnInfo_.names[eq]==fieldName)
       {
           if(eqnInfo_.rank[eq]==SCALAR)
           {
              for(unsigned int id=0; id<scalarFields_.size();id++)
              {
                if(fieldName==scalarFields_[id]->name())
                {
                  solveEquation(scalarFields_[id], eq);
                  return;
                }
              }
           }
           else if(eqnInfo_.rank[eq]==VECTOR)
           {
              for(unsigned int id=0; id<vectorFields_.size();id++)
              {
                if(fieldName==vectorFields_[id]->name())
                {
                  solveEquation(vectorFields_[id], eq);
                  return;
                }
              }
           }
           else if(eqnInfo_.rank[eq]==TENSOR)
           {
              for(unsigned int id=0; id<tensorFields_.size();id++)
              {
                if(fieldName==tensorFields_[id]->name())
                {
                  solveEquation(tensorFields_[id], eq);
                  return;
                }
              }
           }

       }
   }

   FatalError << "Cannot find equation for field " << fieldName;

}

template<class T>
tmp<fvMatrix<T>> AuxEquations::createEquation(
                                GeometricField<T,fvPatchField,volMesh>& field,
                                int eqId
                               )
{
  if(eqnInfo_.type[eqId]==TRANSPORT)
  {
    return transportEquation(field,eqId);
  }

  return equilibriumEquation(field,eqId);
}

template<class rank>
void  AuxEquations::solveEquation(
                         GeometricField<rank,fvPatchField,volMesh> * field,
                         int eqId
                        )
{
  //-Generate base system of equations
  fvMatrix< rank > fieldEqn
  (
    createEquation((*field),eqId)
  );

  //Apply closures
  eqnInfo_.closures[eqId]->closeEquation(fieldEqn,(*field));

  //Solve based on type
  if(eqnInfo_.type[eqId]==TRANSPORT)
  {
    solveDifferentialEquation(fieldEqn,(*field));
  }
  else if(eqnInfo_.type[eqId]==EQUILIBRIUM)
  {
    solveAlgebraicEquation(fieldEqn,(*field));
  }


}

template<class rank>
void AuxEquations::solveDifferentialEquation(
                                             fvMatrix<rank>& diffEqn,
                                             GeometricField<rank,fvPatchField,volMesh>& field
                                            )
{

  diffEqn.relax();
  field.correctBoundaryConditions();
  diffEqn.solve();

}


template<class rank>
void AuxEquations::solveAlgebraicEquation(
                                          fvMatrix<rank>& algEqn,
                                          GeometricField<rank,fvPatchField,volMesh>& field
                                          )
{

 Info << "Solving Equilibrium equation for " << field.name();
 field = algEqn.H()/algEqn.A();

}

const volScalarField* AuxEquations::getFieldScalar(word fieldname)
{
    for(unsigned int id=0;id<scalarFields_.size();id++)
    {
        if(fieldname==scalarFields_[id]->name())
         return  scalarFields_[id];
    }

    FatalError << "Cannot find field " << fieldname;
    return NULL;
}

const volVectorField* AuxEquations::getFieldVector(word fieldname)
{
    for(unsigned int id=0;id<vectorFields_.size();id++)
    {
        if(fieldname==vectorFields_[id]->name())
         return  vectorFields_[id];
    }

    FatalError << "Cannot find field " << fieldname;
    return NULL;
}

const volTensorField* AuxEquations::getFieldTensor(word fieldname)
{
    for(unsigned int id=0;id<tensorFields_.size();id++)
    {
        if(fieldname==tensorFields_[id]->name())
         return  tensorFields_[id];
    }

    FatalError << "Cannot find field " << fieldname;
    return NULL;
}
