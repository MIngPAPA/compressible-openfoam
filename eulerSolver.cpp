/*---------------------------------------------------------------------------*\
    eulerSolver - A simple parallel first order Gas-dynamics solver
                  based on the OpenFOAM library
    Copyright (C) 2012, Pavanakumar Mohanamuraly

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include<fvCFD.H>
#include<UList.H>
#include "fluxSchemes.H"

int main(int argc, char *argv[])
{
   #include "setRootCase.H"
   #include "createTime.H"
   #include "createMesh.H"
   #include "setFreeStream.H"
   #include "createFields.H"
   #include "readFluxScheme.H"

   // Unit face normals
   surfaceVectorField nf = mesh.Sf() / mesh.magSf();
   long int iter = 1;
   scalar rhoResidMax = 0.0;
 
   /// Time step loop
   while( runTime.loop() ) {

     /// Post the non-blocking send/recv of fields
//     #include "initEvaluateFields.H"

     /// Construct the fluxes at faces
     #include "constructFaceFlux.H"
     
     /// Sum-up face fluxes of all associated faces of a cell
     #include "sumFlux.H"
     
     /// Obtain boundary fluxes
     #include "boundaryFlux.H"
     
     /// State update
     #include "stateUpdateLTS.H"

     Info << "Iteration = " << iter++ << " Max residue = " << rhoResidMax << endl;
     /// Solution output
     runTime.write();
   }
 
   return 0;
}


