function maxcalc = Func_Max(y)

syms t

f2 = diff(y,t)==0;
disp(f2);
extreme_points = solve(f2,t);
extreme_values = subs(y, t, extreme_points);
[maxX, maxidx] = max(extreme_values);
maxcalc = simplify(maxX, 'steps', 50);