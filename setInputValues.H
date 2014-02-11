   /// Default values
   scalar M_inf = 0.5, aoa = 0.0, CFL = 0.8;

   /// The inputs are parsed from the command line
   if(  args.optionFound("M") == false ) /// This is a needed input so throw an error
     Info << "Need a Mach number input " << Foam::FatalError.exit();
   else
     M_inf = args.optionRead("M");
   /// Angle of attack input
   if( args.optionFound("aoa") == true )
     aoa = args.optionRead("aoa");
   /// CFL input
   if( args.optionFound("CFL") == true )
     CFL = args.optionRead("CFL");


