/// Default flux solver is roe scheme 
word fluxScheme("roe");
fluxSolver = &roe;

if (mesh.schemesDict().found("fluxScheme"))
{
  fluxScheme = word(mesh.schemesDict().lookup("fluxScheme"));
  if ((fluxScheme == "roe") ){
    Info<< "fluxScheme: " << fluxScheme << endl;
    fluxSolver = &roe;
  }
  else if( (fluxScheme == "van_leer") ){
    Info<< "fluxScheme: " << fluxScheme << endl;
    fluxSolver = &van_leer;
  }
  else{
    FatalErrorIn("EulerSolver::readFluxScheme") << "fluxScheme: " << fluxScheme
                                                << " is not a valid choice. "
                                                << "Options are: roe,van_leer"
                                                << abort(FatalError);
  }
}


