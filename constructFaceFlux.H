
  /// Loop over internal faces and get face flux from L and R cell center
  label leftCell , rightCell;
  forAll( mesh.owner() , iface ) {
    /// Get the left and right cell index
    leftCell = mesh.owner()[iface];
    rightCell = mesh.neighbour()[iface]; 
    /// Approximate Riemann solver at interface 
    (*fluxSolver)( &rho[leftCell]  , &U[leftCell]  , &p[leftCell] , /// Left
         &rho[rightCell] , &U[rightCell] , &p[rightCell] ,          /// Right
         &massFlux[iface] , &momFlux[iface] , &energyFlux[iface] ,
         &nf[iface] ); 
    /// Multiply with face area to get face flux
    massFlux[iface] *= mesh.magSf()[iface];
    momFlux[iface] *= mesh.magSf()[iface];
    energyFlux[iface] *= mesh.magSf()[iface];
  }


