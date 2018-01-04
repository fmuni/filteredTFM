Info<< "Constructing momentum equations" << endl;

MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);
MRF.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

volScalarField KdAve(fluid.KdAve());

{
   // volScalarField Vm(fluid.Vm());

    {
        U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1) + fvm::div(alphaRhoPhi1, U1)
          - fvm::Sp(contErr1, U1)
          + MRF.DDt(alpha1*rho1, U1)
          + phase1.stressClosure().divDevRhoReff(U1)
         ==
           fvOptions(alpha1, rho1, U1)
        );
        U1Eqn.relax();
        fluid.applyDrag(U1,U1Eqn);
        fvOptions.constrain(U1Eqn);
        U1.correctBoundaryConditions();
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2) + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(contErr2, U2)
          + MRF.DDt(alpha2*rho2 , U2)
          + phase2.stressClosure().divDevRhoReff(U2)
         ==
           fvOptions(alpha2, rho2, U2)
        );
        U2Eqn.relax();
        fluid.applyDrag(U2,U2Eqn);
        fvOptions.constrain(U2Eqn);
        U2.correctBoundaryConditions();
        fvOptions.correct(U2);
    }
}
