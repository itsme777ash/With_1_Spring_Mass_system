function fe = Load_vec_matrix(q,x)

le = x(2) - x(1);
fe = [q*le/2;q*le*le/12;q*le/2;-q*le*le/12];
