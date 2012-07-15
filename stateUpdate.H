// Convert to conservative
forAll( rho , icell ) {
  convertConservative( &rho[icell] , &U[icell] , &p[icell] );
}

forAll( rho , icell ) {
  scalar inv_vol = 1.0 / mesh.V()[icell];
  rho[icell] = rho[icell] - inv_vol * runTime.deltaT().value() * massResidue[icell];
  U[icell]   = U[icell]   - inv_vol * runTime.deltaT().value() * momResidue[icell];
  p[icell]   = p[icell]   - inv_vol * runTime.deltaT().value() * energyResidue[icell];
}

// Convert to primitive
forAll( rho , icell ) {
  convertPrimitive( &rho[icell] , &U[icell] , &p[icell] );  
}
 
