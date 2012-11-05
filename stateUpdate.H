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

  /// Convert to conservative
  forAll( rho , icell )
    convertConservative( &rho[icell] , &U[icell] , &p[icell] );

  /// Euler state update
  forAll( rho , icell ) {
    scalar inv_vol = 1.0 / mesh.V()[icell];
    rho[icell] = rho[icell] - inv_vol * runTime.deltaT().value() * massResidue[icell];
    U[icell]   = U[icell]   - inv_vol * runTime.deltaT().value() * momResidue[icell];
    p[icell]   = p[icell]   - inv_vol * runTime.deltaT().value() * energyResidue[icell];
  }

  /// Convert to primitive
  forAll( rho , icell )
    convertPrimitive( &rho[icell] , &U[icell] , &p[icell] );
 
