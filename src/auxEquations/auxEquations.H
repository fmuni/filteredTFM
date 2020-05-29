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

#ifndef AuxEquations_H
#define AuxEquations_H

#include "fvCFD.H"
#include "twoPhaseSystem.H"


namespace Foam
{

enum
{
    SCALAR, VECTOR, TENSOR
};

enum
{
    TRANSPORT, EQUILIBRIUM
};

enum
{
    MIXTURE, CONTINUOUS, DISPERSED
};

class EqnClosure;

//-Keep generic information for easy addressing
struct eqnInformation
{
    wordList                    names; //names of fields
    labelList                    rank; //rank (SCALAR,VECTOR or TENSOR)
    labelList                    type; //transport or algebraic
    PtrList<EqnClosure>      closures; //Closure for the correspondinf equation
    labelList                   phase; //phase knowledge is required for transport
};




class AuxEquations
{
    protected:

    //Members

    twoPhaseSystem&               fluid_;

    dictionary                 settings_;

    eqnInformation              eqnInfo_;

    //Field containers
    PtrList<volScalarField>   scalarFields_;
    PtrList<volVectorField>   vectorFields_;
    PtrList<volTensorField>   tensorFields_;

    //Methods

    //-Creates a transport equation
    template<class T>
    tmp<fvMatrix<T>> transportEquation
    (
        GeometricField<T,fvPatchField,volMesh>& field,
        label eqId
    );

    //-Creates an equilibrium equation
    template<class T>
    tmp<fvMatrix<T>> equilibriumEquation
    (
        GeometricField<T,fvPatchField,volMesh>& field,
        int eqId
    );

    //-Creates system of equations
    template<class T>
    tmp<fvMatrix<T>> createEquation
    (
        GeometricField<T,fvPatchField,volMesh>& field,
        int eqId
    );

    //-Build and solve equations
    template<class rank>
    void solveEquation
    (
        GeometricField<rank,fvPatchField,volMesh>& field,
        int eqId
    );

    template<class rank>
    void solveDifferentialEquation
    (
        fvMatrix<rank>& diffEqn,
        GeometricField<rank,fvPatchField,volMesh>& field
    );

    template<class rank>
    void solveAlgebraicEquation
    (
        fvMatrix<rank>& algEqn,
        GeometricField<rank,fvPatchField,volMesh>& field
    );


    public:


    //Members



    //Methods

    //Construct from fTFM object and dictionary
    AuxEquations
    (
        twoPhaseSystem&             fluid,
        dictionary              auxEquationsDict
    );

    //Destructor
    ~AuxEquations();

    void requestEquation(word equationName, word requestingClass) const;

    //Check field and call equation builder and solver
    void solveField(word fieldName);

    //Getters for fieldsToSolve
    const volScalarField& getFieldScalar(word fieldname);
    const volVectorField& getFieldVector(word fieldname);
    const volTensorField& getFieldTensor(word fieldname);

    //Getters for id
    label getEquationId(word fieldName);

};
}

#endif
