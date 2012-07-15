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



