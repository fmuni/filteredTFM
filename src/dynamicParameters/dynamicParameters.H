#ifndef DYNAMIC_PARAMETERS_H
#define DYNAMIC_PARAMETERS_H
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
#include "fvCFD.H"
#include "multiphaseFlowBasic.H"


namespace Foam
{
 //-Forward class declarations
 class twoPhaseSystem;
 class phaseModel;

 enum {DIMENSIONLESS, DIMENSIONAL};

 class DynamicParameters
 {
  protected:

   //Members


    const twoPhaseSystem&             fluid_;

    dictionary                     settings_;

    multiphaseFlowBasic   basicCalculations_;

    //-Reference filter size
    scalar      refDeltaf_;

    //-Reference maximum particle volume fraction
    scalar        phiPMax_;

    //-filter width
    scalar        filterWidth_;

    //-Important markers

     //-scaled slip Velocity
     mutable volVectorField slipU_;

     //-strain rate
     mutable volScalarField Sr_;


   //Methods

   //-Calculate quantities for an isolated sphere
   void setupSettling(word dragLaw, double filterSize);

   //-Check if the mesh is Cartesian
   void checkCartesian();



   //-Returns filter Freud number
   double  (DynamicParameters::*FrFilter)(int celli) const;

   //-Returns filter size
   double  (DynamicParameters::*deltaF)(int celli, int dim) const;

   //Return Filter Freud number Cartesian
   double FrConstant(int celli) const;

   //Return filter Freud number non-Cartesian
   double FrNoConstant(int celli) const;

   //Return filter size Cartesian
   double deltaFConstant(int celli, int dim) const;

   //Return filter size non-Cartesian
   double deltaFNoConstant(int celli, int dim) const;


   //Vector of processed scalar fields
   //i.e., some operation has been performed on them
   std::vector<volScalarField*>    filteredFieldsScalar_;

  public:


    //Members

    //Methods

     //Construct from fTFM object and dictionary
     DynamicParameters
     (
      const twoPhaseSystem&             fluid,
      dictionary              auxEquationsDict
     );

     //Destructor
     ~DynamicParameters();

     //Return filter Freud number
     double filterFr(int celli) const;

     //Get desired field
     //Base-field + operation name
     const volScalarField& getProcessedField(
                                              volScalarField& baseField,
                                              word operationName
    ) const;

    //Field name + operation name
    //Base-field + operation name
    const volScalarField& getProcessedField(
                                             word fieldName,
                                             word operationName
    ) const;

    //Create new operation for scalar field
    void createOperationField(
                                     volScalarField& baseField,
                                     word operationName
    ) const;

    //-Get settling velocity
    double settlingU() const;

    //-Get filter size
    double filterSize( label celli, int dim=DIMENSIONLESS) const;

    //-Get slip velocity
    const volVectorField& scaledSlipU(const phaseModel& continuous,
                                    const  phaseModel& dispersed) const;

    //-Get strain rate
    const volScalarField& strainRate(const phaseModel& phase) const;

    //-Return gravity value
    double gValue() const;

    //- Return gravity vector
    dimensionedVector g()  const;

    //- Return error if reference filter size is larger than actual filter size
    void checkFilterSize(scalar refFilterSize) const;



 };
}

#endif
