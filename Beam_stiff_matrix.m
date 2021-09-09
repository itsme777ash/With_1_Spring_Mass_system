function ke = Beam_stiff_matrix( E, Ie, x)

le = x(2) - x(1);

factor = E*Ie/(le^3);

ke = factor * [12 , 6*le, -12, 6*le;
               6*le, 4*le*le, -6*le, 2*le*le;
               -12, -6*le, 12, -6*le;
               6*le, 2*le*le, -6*le, 4*le*le];
           