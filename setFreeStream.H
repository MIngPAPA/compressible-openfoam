   /// Input M_inf
   scalar M_inf = 2.0, CFL = 0.5;
   /// Freestream values
   scalar rho_inf = 1.0, p_inf = 1.0 / ( gama * M_inf );
   vector u_inf ( 1.0 ,0.0 ,0.0 );
   /// Set the correct values
   rho = rho_inf;
   U   = u_inf;
   p   = p_inf;
 
