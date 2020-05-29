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

using namespace Foam;
//-------------------------- Constructors ---------------------------------//
AuxEquations::AuxEquations
(
    twoPhaseSystem&             fluid,
    dictionary       auxEquationsDict
)
:
fluid_(fluid),
settings_(auxEquationsDict)
{
    //Fill information container
    if(settings_.found("fieldsToSolve"))
    {
        eqnInfo_.names = wordList(settings_.lookup("fieldsToSolve"));
    }

    forAll(eqnInfo_.names,eq)
    {
        dictionary eqnDict(settings_.subDict(eqnInfo_.names[eq]));

        word type(eqnDict.lookup("equation"));

        std::size_t found(0);

        //Check field rank and call sield constructor
        if(type.find("scalar")!=std::string::npos)
        {
            eqnInfo_.rank.append(SCALAR);
            found = type.find("scalar");
            scalarFields_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        eqnInfo_.names[eq],
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluid_.mesh()
                )
            );

        }
        else if(type.find("vector")!=std::string::npos)
        {
            eqnInfo_.rank.append(VECTOR);
            found = type.find("vector");

            vectorFields_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        eqnInfo_.names[eq],
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluid_.mesh()
                )
            );

        }
        else if(type.find("tensor")!=std::string::npos)
        {
            eqnInfo_.rank.append(TENSOR);
            found = type.find("tensor");

            tensorFields_.append
            (
                new volTensorField
                (
                    IOobject
                    (
                        eqnInfo_.names[eq],
                        fluid_.mesh().time().timeName(),
                        fluid_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluid_.mesh()
                )
            );

        }
        else
        {
            FatalErrorInFunction<< "Invalid equation " << type
            << "for field" << eqnInfo_.names[eq]
            << abort(FatalError);
        }

        //Check equation type
        if(type.find("Transport",found+1)!=std::string::npos)
        {
            eqnInfo_.type.append(TRANSPORT);

            //Assign transporting phase
            word trnsPhase(eqnDict.lookup("carrierPhase"));
            if(trnsPhase=="mixture")
            {
                eqnInfo_.phase.append(MIXTURE);
            }
            else if(trnsPhase=="continuous")
            {
                eqnInfo_.phase.append(CONTINUOUS);
            }
            else if(trnsPhase=="dispersed")
            {
                eqnInfo_.phase.append(DISPERSED);
            }
            else
            {
                FatalErrorInFunction<< "Invalid equation " << type
                << "for field" << eqnInfo_.names[eq] << "\n"
                << "Carrier phase is not specified!"
                << abort(FatalError);
            }

        }
        else if(type.find("Equilibrium",found+1)!=std::string::npos)
        {
            eqnInfo_.type.append(EQUILIBRIUM);
            eqnInfo_.phase.append(MIXTURE);
        }
        else
        {
            FatalErrorInFunction<< "Invalid equation " << type
            << "for field" << eqnInfo_.names[eq]
            << abort(FatalError);
        }

        Info << "Found " << type << " equation for field "
        << eqnInfo_.names[eq];

        //Read closure
        eqnInfo_.closures.append
        (
            new EqnClosure
            (
                fluid_,
                settings_,
                eqnInfo_.names[eq]
            )
        );

    }


};
//-------------------------- Destructors ----------------------------------//
AuxEquations::~AuxEquations()
{
};
//---------------------------   Methods  ----------------------------------//
void AuxEquations::requestEquation(word equationName, word requestingClass) const
{
    forAll(eqnInfo_.names,eq)
    {
        if(equationName == eqnInfo_.names[eq])
        return;
    }

    FatalErrorInFunction<< "Cannot find equation " << equationName
                        << "requested by " << requestingClass << " !! "
                        << abort(FatalError);
}

template<class T>
tmp<fvMatrix<T>> AuxEquations::transportEquation
(
    GeometricField<T,fvPatchField,volMesh>& field,
    label eqId
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

    phaseModel* alpha(NULL);

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
      - fvc::Sp
        (
            fvc::ddt((*alpha), alpha->rho())
            + fvc::div(alpha->alphaRhoPhi()),
            field
        )

    );

}

template<class T>
tmp<fvMatrix<T>> AuxEquations::equilibriumEquation
(
    GeometricField<T,fvPatchField,volMesh>& field,
    label eqId
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
    label eq(getEquationId(fieldName));

    if(eqnInfo_.rank[eq]==SCALAR)
    {
        forAll(scalarFields_,id)
        {
            if(fieldName==scalarFields_[id].name())
            {
                solveEquation(scalarFields_[id], eq);
                return;
            }
        }
    }
    else if(eqnInfo_.rank[eq]==VECTOR)
    {
        forAll(vectorFields_,id)
        {
            if(fieldName==vectorFields_[id].name())
            {
                solveEquation(vectorFields_[id], eq);
                return;
            }
        }
    }
    else if(eqnInfo_.rank[eq]==TENSOR)
    {
        forAll(tensorFields_,id)
        {
            if(fieldName==tensorFields_[id].name())
            {
                solveEquation(tensorFields_[id], eq);
                return;
            }
        }
    }

    FatalErrorInFunction<< "Cannot find equation for field " << fieldName
                        << abort(FatalError);

}

template<class T>
tmp<fvMatrix<T>> AuxEquations::createEquation
(
    GeometricField<T,fvPatchField,volMesh>& field,
    label eqId
)
{
    if(eqnInfo_.type[eqId]==TRANSPORT)
    {
        return transportEquation(field,eqId);
    }

    return equilibriumEquation(field,eqId);
}

template<class rank>
void  AuxEquations::solveEquation
(
    GeometricField<rank,fvPatchField,volMesh>& field,
    label eqId
)
{
    //-Generate base system of equations
    fvMatrix<rank> fieldEqn
    (
        createEquation(field,eqId)
    );

    //Apply closures
    eqnInfo_.closures[eqId].closeEquation(fieldEqn,field);

    //Solve based on type
    if(eqnInfo_.type[eqId]==TRANSPORT)
    {
        solveDifferentialEquation(fieldEqn,field);
    }
    else if(eqnInfo_.type[eqId]==EQUILIBRIUM)
    {
        solveAlgebraicEquation(fieldEqn,field);
    }


}

template<class rank>
void AuxEquations::solveDifferentialEquation
(
    fvMatrix<rank>& diffEqn,
    GeometricField<rank,fvPatchField,volMesh>& field
)
{

    diffEqn.relax();
    field.correctBoundaryConditions();
    diffEqn.solve();

}


template<class rank>
void AuxEquations::solveAlgebraicEquation
(
    fvMatrix<rank>& algEqn,
    GeometricField<rank,fvPatchField,volMesh>& field
)
{

    Info << "Solving Equilibrium equation for " << field.name();
    field = algEqn.H()/algEqn.A();

}

const volScalarField& AuxEquations::getFieldScalar(word fieldname)
{
    forAll(scalarFields_,id)
    {
        if(fieldname==scalarFields_[id].name())
        {
            return  scalarFields_[id];
        }
    }

    FatalErrorInFunction<< "Cannot find field " << fieldname
                        << abort(FatalError);
    return scalarFields_[0];
}

const volVectorField& AuxEquations::getFieldVector(word fieldname)
{
    forAll(vectorFields_,id)
    {
        if(fieldname==vectorFields_[id].name())
        {
            return  vectorFields_[id];
        }
    }

    FatalErrorInFunction<< "Cannot find field " << fieldname
                        << abort(FatalError);
    return vectorFields_[0];
}


const volTensorField& AuxEquations::getFieldTensor(word fieldname)
{
    forAll(tensorFields_,id)
    {
        if(fieldname==tensorFields_[id].name())
        {
            return  tensorFields_[id];
        }
    }

    FatalErrorInFunction<< "Cannot find field " << fieldname
                        << abort(FatalError);
    return tensorFields_[0];
}

label  AuxEquations::getEquationId(word fieldname)
{
    forAll(eqnInfo_.names,eq)
    {
        if(eqnInfo_.names[eq]==fieldname)
        {
            return eq;
        }
    }

    FatalErrorInFunction<< "Cannot find equation " << fieldname
                        << abort(FatalError);
    return 0;
}
