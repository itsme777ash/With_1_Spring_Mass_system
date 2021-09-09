function Mass = Beam_Mass_matrix(Area , density , l)

% Mass = density*Area*l*0.5*[1,0,0,0;
%                            0,0,0,0;
%                            0,0,1,0;
%                            0,0,0,0];
factor = density*Area*l/420;
Mass = factor*[156,22*l,54,-13*l;
               22*l,4*l*l,13*l,-3*l*l;
               54,13*l,156,-22*l;
               -13*l,-3*l*l,-22*l,4*l*l];