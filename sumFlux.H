
   /// Clear all residues
   massResidue = dimensionedScalar( "", massResidue.dimensions() , 0.0 );
   momResidue = dimensionedVector( "" , momResidue.dimensions() , vector( 0.0 , 0.0 , 0.0 ) );
   energyResidue = dimensionedScalar( "" , energyResidue.dimensions() ,  0.0 );
   /// Loop over each face and add the flux to left and subtract from the right
   forAll( mesh.owner() , iface ) {
     /// Store left and right cell reference
     leftCell = mesh.owner()[iface];
     rightCell = mesh.neighbour()[iface];
     /// Note that the normal vector to a face will point
     /// from the owner cell to the neighbour cell (L to R)
     /// Add to left cell
     massResidue[leftCell]    += massFlux[iface];
     momResidue[leftCell]     += momFlux[iface];
     energyResidue[leftCell]  += energyFlux[iface];
     /// Subtract from right cell
     massResidue[rightCell]   -= massFlux[iface];
     momResidue[rightCell]    -= momFlux[iface];
     energyResidue[rightCell] -= energyFlux[iface];
   }
