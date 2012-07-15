/*---------------------------------------------------------------------------*\
Author
   Pavanakumar Mohanamuraly
   Scientist Fellow, CTFD division,
   NAL, CSIR

Application
    euler_solver

Description
   A simple parallel first order Gas-dynamics solver

WARNING:
   This material is for use in OpenFOAM discussions held at NAL CTFD division. 
   Not for distribution outside of CSIR.  
\*---------------------------------------------------------------------------*/

#include<fvCFD.H>
#include<UList.H>

/// Riemann solver that takes in density,velocity and 
/// pressure and return the flux
void (*fluxSolver)( scalar *rhoL , vector *uL , scalar *pL , 
          scalar *rhoR , vector *uR , scalar *pR ,
          scalar *mass , vector *mom , scalar *energy ,
          vector *normal ) = NULL;
/// Roe Flux scheme
void roe( scalar *rhoL , vector *uL , scalar *pL , 
          scalar *rhoR , vector *uR , scalar *pR ,
          scalar *mass , vector *mom , scalar *energy ,
          vector *normal );


/// Returns the normal flux given the state vector
void normalFlux( scalar *rho , vector *u , scalar *p , 
                 vector *normal , scalar flux[5] );

void convertConservative( scalar *rho , vector *u , scalar *p );
void convertPrimitive( scalar *mass , vector *mom , scalar *energy );
 
/// Constants
scalar gama = 1.4;
scalar gbygm1 = gama / ( gama - 1.0 ) , gm1 = gama - 1.0 ;

int main(int argc, char *argv[])
{
   #include "setRootCase.H"
   #include "createTime.H"
   #include "createMesh.H"
   #include "createFields.H"

   /// Freestream values
   scalar rho_inf = 1.228 , p_inf = 101325.0, M_inf = 2.0 ;
   scalar magU_inf = std::sqrt( 1.40e0 * p_inf / rho_inf ) * M_inf;
   vector u_inf ( magU_inf , 0.0 , 0.0 ) ;
   // Unit face normals
   surfaceVectorField nf = mesh.Sf() / mesh.magSf();
   //nf.boundaryField() = mesh.Sf().boundaryField() / mesh.magSf().boundaryField();
   label iter = 1;
   fluxSolver = &roe;  /// Use roe flux solver

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

/// Try implementing your own riemann flux scheme
/// This is Roe's Approximate Riemann solver
void roe( scalar *rhoL , vector *uL , scalar *pL ,
          scalar *rhoR , vector *uR , scalar *pR ,
          scalar *mass , vector *mom , scalar *energy ,
          vector *normal ){
  scalar h0L, h0R, Rfac, rhoRoe, aRoe,
         h0Roe, drho, dp, VnL, VnR, VnRoe, qSqRoe,
         dV, drhoTerm , uL2 , uR2;
  scalar fL[5], fR[5], df1[5], df2[5], df3[5] , dot_prod;
  vector dvel , velRoe;
  //  Left and Right Flow values
  uL2 = (*uL) & (*uL);
  uR2 = (*uR) & (*uR);
  h0L = gbygm1 * (*pL) / (*rhoL) + 0.50 * ( uL2 );
  h0R = gbygm1 * (*pR) / (*rhoR) + 0.50 * ( uR2 );
  VnL = (*uL) & (*normal);
  VnR = (*uR) & (*normal);
  //  Roe Averaged States
  Rfac = std::sqrt( (*rhoR) / (*rhoL) );
  rhoRoe = Rfac * (*rhoL);
  velRoe = ( (*uL) + Rfac * (*uR) ) / ( 1.0 + Rfac );
  h0Roe = (h0L + Rfac * h0R) / ( 1.0 + Rfac );
  qSqRoe = velRoe & velRoe ;
  aRoe = std::sqrt( gm1 * ( h0Roe - 0.50 * qSqRoe ) );
  VnRoe = (VnL + Rfac * VnR) / (1.0 + Rfac);
  //  Differential values for characteristics
  drho = (*rhoR) - (*rhoL);
  dvel = (*uR) - (*uL);
  dp   = (*pR) - (*pL);
  dV = dvel & (*normal);
  normalFlux( rhoL , uL , pL , normal , fL );
  normalFlux( rhoR , uR , pR , normal , fR );
  // Form Inteface roe averaged fluxes df1 to df3
  // Calculate df1
  drhoTerm = drho - dp / ( aRoe * aRoe );
  df1[0]   = std::abs(VnRoe) * drhoTerm;
  dot_prod = 0.0;
  for( int i = 0 ; i < 3 ; ++i ){
    df1[i+1] = std::abs(VnRoe) * ( velRoe[i] * drhoTerm +
               rhoRoe * (dvel[i] - (*normal)[i] * dV) );
    dot_prod += velRoe[i] * dvel[i];
  }
  df1[4]   = std::abs(VnRoe) * ( 0.5 * qSqRoe * drhoTerm +
             rhoRoe * ( dot_prod - VnRoe * dV ) );
  // Calculate df2
  df2[0]   = 0.5 * std::abs(VnRoe + aRoe) *
             ( dp / ( aRoe * aRoe ) + rhoRoe * dV / aRoe );
  for( int i = 0 ; i < 3 ; ++i )
    df2[i+1] = df2[0] * ( velRoe[i] + (*normal)[i] * aRoe );
  df2[4] = df2[0] * ( h0Roe + VnRoe * aRoe );
  // Calculate df3
  df3[0]   = 0.5 * std::abs(VnRoe - aRoe) *
             ( dp / ( aRoe * aRoe ) - rhoRoe * dV / aRoe );
  for( int i = 0 ; i < 3 ; ++i )
    df3[i+1] = df3[0] * ( velRoe[i] - (*normal)[i] * aRoe );
  df3[4]   = df3[0] * ( h0Roe - VnRoe * aRoe );
  // Sum fluxes 0.5 * (fL + fR - df1 - df2 -df3 ) to yield flux at the interface
  (*mass) = 0.5 * ( fL[0] + fR[0] - df1[0] - df2[0] - df3[0] );
  for( int i = 1 ; i < 4 ; ++i )
    (*mom)[i-1] = 0.5 * ( fL[i] + fR[i] - df1[i] - df2[i] - df3[i] );
  (*energy) = 0.5 * ( fL[4] + fR[4] - df1[4] - df2[4] - df3[4] );
}

/// The normal flux to a face given the state vector 
/// and normal at a face
void normalFlux( scalar *rho , vector *u , scalar *p ,
                 vector *normal , scalar f[5] ){
  f[0] = (*rho) * ( (*u) & (*normal) );
  for( int i = 0 ; i < 3 ; ++i )
    f[i+1] = f[0] * (*u)[i] + (*p) * (*normal)[i];
  f[4] = f[0] * ( gbygm1 * (*p) / (*rho) + 0.5 * ( (*u) & (*u) ) );
}

//
void convertConservative( scalar *rho , vector *u , scalar *p ){
  // \rho e_t = p / ( \gamma - 1 ) + 1/2 * \rho * u^2
  (*p) = (*p) / gm1 + 0.5 * (*rho) * ( (*u) & (*u) );
  // (\rho u,\rho v,\rho w) = (u,v,w) * \rho
  (*u) *= (*rho);
}

//
void convertPrimitive( scalar *mass , vector *mom , scalar *energy ){
  (*mom) /= (*mass);
  // p = (\gamma - 1) (\rho e_t - 1/2 \rho u^2 )
  (*energy) = gm1 * ( (*energy) - 0.5 * (*mass) * ( (*mom) & (*mom) ) );
}

