   /// Remember BC's for Gas Dynamics solvers
   /// are based on characteristics and cannot be applied
   /// to individual variables. Rather it should be applied 
   /// to the state vector as a whole. OpenFOAM framework 
   /// does not allow us to do this without breaking the 
   /// convention.
   
   /// Loop over boundary patches and determine the 
   /// type of boundary condition of the patch 
   forAll ( mesh.boundaryMesh() , ipatch ) {
     scalar bflux[5];
     word BCtype = mesh.boundaryMesh().physicalTypes()[ipatch];  
     word BCname = mesh.boundaryMesh().names()[ipatch];
     const UList<label> &bfaceCells = mesh.boundaryMesh()[ipatch].faceCells();
/////////////////////////////////////////////////////////////////////////////////
                                /// Slip Wall BC ///
/////////////////////////////////////////////////////////////////////////////////
     if( BCtype == "slip" ) {
       forAll( bfaceCells  , iface ) {
         /// Extrapolate wall pressure
         scalar p_e = p[ bfaceCells[iface] ]  ;
         vector normal = nf.boundaryField()[ipatch][iface]; 
         scalar face_area = mesh.magSf().boundaryField()[ipatch][iface];
         momResidue[ bfaceCells[iface] ] += p_e * normal * face_area;
       }
     }
/////////////////////////////////////////////////////////////////////////////////
                          /// Extrapolation Outlet BC /// 
/////////////////////////////////////////////////////////////////////////////////
     if( BCtype == "extrapolatedOutflow" ) {
       forAll( bfaceCells  , iface ) {
         label myCell = bfaceCells[iface];
         vector normal = nf.boundaryField()[ipatch][iface]; 
         scalar face_area = mesh.magSf().boundaryField()[ipatch][iface];
         /// Use adjacent cell center value of face to get flux (zeroth order)
         normalFlux( &rho[myCell] , &U[myCell] , &p[myCell] , &normal , bflux );
         massResidue[ myCell ]   += bflux[0] * face_area;
         momResidue[ myCell ][0] += bflux[1] * face_area;
         momResidue[ myCell ][1] += bflux[2] * face_area;
         momResidue[ myCell ][2] += bflux[3] * face_area;
         energyResidue[ myCell ] += bflux[4] * face_area;
       }
     }
/////////////////////////////////////////////////////////////////////////////////
                     /// Supersonic Inlet BC ( Drichlet BC )  ///
/////////////////////////////////////////////////////////////////////////////////
     if( BCtype == "supersonicInlet" ) {
       forAll( bfaceCells  , iface ) {
         label myCell = bfaceCells[iface];
         vector normal = nf.boundaryField()[ipatch][iface] ;
         scalar face_area = mesh.magSf().boundaryField()[ipatch][iface];
         /// Use free-stream values to calculate face flux
         normalFlux( &rho_inf , &u_inf , &p_inf , &normal , bflux );
         massResidue[ myCell ]   += bflux[0] * face_area;
         momResidue[ myCell ][0] += bflux[1] * face_area;
         momResidue[ myCell ][1] += bflux[2] * face_area;
         momResidue[ myCell ][2] += bflux[3] * face_area;
         energyResidue[ myCell ] += bflux[4] * face_area;
       }       
     }     
     /// Processor BC (parallel)
     
   }

