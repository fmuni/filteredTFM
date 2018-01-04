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
#include <vector>


namespace Foam
{

 enum{
       SCALAR, VECTOR, TENSOR
     };

 enum{
       TRANSPORT, EQUILIBRIUM
     };

 enum{
       MIXTURE, CONTINUOUS, DISPERSED
     };

  class EqnClosure;

  //-Keep generic information for easy addressing
  struct eqnInformation
  {
     wordList                    names; //names of fields
     std::vector<int>             rank; //rank (SCALAR,VECTOR or TENSOR)
     std::vector<int>             type; //transport or algebraic
     std::vector<EqnClosure*> closures; //Closure for the correspondinf equation
     std::vector<int>            phase; //phase knowledge is required for transport
  };




 class AuxEquations
 {
  protected:

   //Members

    twoPhaseSystem&         fluid_;

    dictionary          settings_;

    eqnInformation              eqnInfo_;


    //Field containers
    std::vector<volScalarField*>   scalarFields_;
    std::vector<volVectorField*>   vectorFields_;
    std::vector<volTensorField*>   tensorFields_;

   //Methods

   //-Creates a transport equation
   template<class T>
   tmp<fvMatrix<T>> transportEquation(
                                       GeometricField<T,fvPatchField,volMesh>& field,
                                       int eqId
                                      );

   //-Creates an equilibrium equation
   template<class T>
   tmp<fvMatrix<T>> equilibriumEquation(
                                         GeometricField<T,fvPatchField,volMesh>& field,
                                         int eqId
                                        );
    //-Creates system of equations
   template<class T>
   tmp<fvMatrix<T>> createEquation(
                                     GeometricField<T,fvPatchField,volMesh>& field,
                                     int eqId
                                   );



    //-Builds and solves equations
    template<class rank>
    void solveEquation(
                             GeometricField<rank,fvPatchField,volMesh> * field,
                             int eqId
                           );

    template<class rank>
    void solveDifferentialEquation(
                                    fvMatrix<rank>& diffEqn,
                                    GeometricField<rank,fvPatchField,volMesh>& field
                                  );

    template<class rank>
    void solveAlgebraicEquation(
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
     const volScalarField* getFieldScalar(word fieldname);
     const volVectorField* getFieldVector(word fieldname);
     const volTensorField* getFieldTensor(word fieldname);

 };
}

#endif
