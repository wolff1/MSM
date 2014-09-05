reset
set term aqua enhanced
set termoption dashed
set xlabel 'u'
set ylabel 'f(u)'
set style data lines
set yrange [   0.000 :   0.833 ]
plot data_file using 1:2 with lines title "{/Symbol F}(2u-2)", "" using 1:3 with lines title "{/Symbol F}(2u-1)", "" using 1:4 with lines title "{/Symbol F}(2u+0)", "" using 1:5 with lines title "{/Symbol F}(2u+1)", "" using 1:6 with lines title "{/Symbol F}(2u+2)", "" using 1:7 with lines title "{/Symbol S}{/Symbol F}(2u+v)" lc rgb "red", "" using 1:8 with lines title "{/Symbol F}(u)" lc rgb "black"
pause -1
quit
