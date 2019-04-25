clear;set terminal qt font 'Helvetica, 24';set parametric;
set isosample 100;
set vrange [-10:10];
set xrange [-10:10];set yrange [-10:10];set zrange [-10:10];
sp 1.000000 + 0.000000*cos(u) + -2.000000*sin(u) + 0.000000*(abs(v) < 5.000000 ? v:NaN), -4.000000 + 2.000000*cos(u) + 0.000000*sin(u) + 0.000000*(abs(v) < 5.000000 ? v:NaN), 0.000000 + 0.000000*cos(u) + 0.000000*sin(u) + 1.000000*(abs(v) < 5.000000 ? v:NaN) title 'P1', 0.000000 + 0.000000*cos(u) + -1.490712*sin(u) + 0.666667*(abs(v) < 5.000000 ? v:NaN), -10.000000 + 0.894427*cos(u) + 1.192570*sin(u) + 0.666667*(abs(v) < 5.000000 ? v:NaN), 2.000000 + -1.788854*cos(u) + 0.596285*sin(u) + 0.333333*(abs(v) < 5.000000 ? v:NaN) title 'P2' ;set title 'overlap 3'
pause -1
