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

//////////////////////
/// Fluid Density  ///
//////////////////////
volScalarField rho
(
  IOobject
  (
    "rho",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);
//////////////////////
/// Fluid Velocity ///
//////////////////////
volVectorField U
(
  IOobject
  (
    "U",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);
//////////////////////
/// Fluid Pressure ///
//////////////////////
volScalarField p
(
  IOobject
  (
    "p",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);
//////////////////////////////////////////////
/// Face Momentum Flux ( \rhouu \cdot ds ) ///
//////////////////////////////////////////////
surfaceVectorField momFlux
(
    IOobject
    (
        "momFlux",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector( "momFlux", 
       rho.dimensions() * U.dimensions() * U.dimensions() * mesh.Sf().dimensions(), 
       vector( 0.0 , 0.0 , 0.0 ) )
);
//////////////////////////////////////////
/// Face Mass Flux ( \rho u \cdot ds ) ///
//////////////////////////////////////////
surfaceScalarField massFlux
(
    IOobject
    (
        "massFlux",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "massFlux", rho.dimensions() * U.dimensions() * mesh.Sf().dimensions(), 0.0 )
);
////////////////////////
/// Face Energy Flux ///
////////////////////////
surfaceScalarField energyFlux
(
    IOobject
    (
        "energyFlux",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "energyFlux", p.dimensions() * U.dimensions() * mesh.Sf().dimensions(), 0.0 )
);
//////////////////////
///  Mass residue  ///
//////////////////////
volScalarField massResidue
(
    IOobject
    (
        "massResidue",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "massResidue", 
       rho.dimensions() * U.dimensions() * U.dimensions() * mesh.Sf().dimensions(), 0.0)
);
////////////////////////
/// Momentum residue ///
////////////////////////
volVectorField momResidue
(
    IOobject
    (
        "momResidue",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector( "momResidue", 
       rho.dimensions() * U.dimensions() * U.dimensions() * mesh.Sf().dimensions(), 
       vector( 0.0 , 0.0 , 0.0 ) )
);
//////////////////////
/// Energy residue ///
//////////////////////
volScalarField energyResidue
(
    IOobject
    (
        "energyResidue",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "energyResidue", p.dimensions() * U.dimensions() * mesh.Sf().dimensions(), 0.0 )
);
///////////////////////
/// Local Time Step ///
///////////////////////
volScalarField localDt
(
    IOobject
    (
        "localDt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar( "localDt", dimless , 0.0 )
);

