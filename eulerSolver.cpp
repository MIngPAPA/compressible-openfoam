/*---------------------------------------------------------------------------*\
Author
   Pavanakumar Mohanamuraly
   Scientist Fellow, CTFD division,
   NAL, CSIR

Application
    eulerSolver

Description
   A simple parallel first order Gas-dynamics solver

WARNING:
   This material is for use in OpenFOAM discussions held at NAL CTFD division. 
   Not for distribution outside of CSIR.
\*---------------------------------------------------------------------------*/

#include<fvCFD.H>
#include<UList.H>
#include "fluxSchemes.H"

int main(int argc, char *argv[])
{
   #include "setRootCase.H"
   #include "createTime.H"
   #include "createMesh.H"
   #include "createFields.H"
   #include "readFluxScheme.H"

   /// Freestream values
   scalar rho_inf = 1.228 , p_inf = 101325.0, M_inf = 2.0 ;
   scalar magU_inf = std::sqrt( 1.40e0 * p_inf / rho_inf ) * M_inf;
   vector u_inf ( magU_inf , 0.0 , 0.0 );
   // Unit face normals
   surfaceVectorField nf = mesh.Sf() / mesh.magSf();
   label iter = 1;

   /// Time step loop
   while( runTime.loop() ) {
     Info<< "Iteration = " << iter++ << nl << endl;

     /// Construct the fluxes at faces
     #include "constructFaceFlux.H"
     
     /// Sum-up face fluxes of all associated faces of a cell
     #include "sumFlux.H"
     
     /// Obtain boundary fluxes
     #include "boundaryFlux.H"
     
     /// State update
     #include "stateUpdate.H"

     // Solution output
     runTime.write();
   }
 
   return 0;
}


