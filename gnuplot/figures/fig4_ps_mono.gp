reset
#set term aqua
#set termoption dash
#set term epslatex mono
set term postscript eps enhanced mono
set output "bwfig4.eps"
#set xlabel '$u$'
#set ylabel '$f(u)$'
set xlabel 'u'
set ylabel 'f(u)'
set style data lines
set yrange [   0.000 :   0.833 ]
#set key at 0.10,0.50 bottom left nobox
#plot data_file using 1:2 with lines title "Phi(2u-2)" ls 10, "" using 1:3 with lines title "Phi(2u-1)" ls 15, "" using 1:4 with lines title "Phi(2u+0)" ls 4, "" using 1:5 with lines title "Phi(2u+1)" ls 24, "" using 1:6 with lines title "Phi(2u+2)" ls 19, "" using 1:7 with lines title "Sum of fine grids" ls 12, "" using 1:8 with lines title "Phi(u)" ls 16
#plot data_file using 1:2 with lines title "Phi(2u-2)" ls 7, "" using 1:3 with lines title "Phi(2u-1)" ls 9, "" using 1:4 with lines title "Phi(2u+0)" ls 5, "" using 1:5 with lines title "Phi(2u+1)" ls 2, "" using 1:6 with lines title "Phi(2u+2)" ls 3, "" using 1:7 with lines title "Sum of fine grids" ls 8, "" using 1:8 with lines title "Phi(u)" ls 4
#plot data_file using 1:2 with lines title "Phi(2u-2)", "" using 1:3 with lines title "Phi(2u-1)", "" using 1:4 with lines title "Phi(2u+0)", "" using 1:5 with lines title "Phi(2u+1)", "" using 1:6 with lines title "Phi(2u+2)", "" using 1:7 with lines title "Sum of fine grids", "" using 1:8 with lines title "Phi(u)"
#plot data_file using 1:2 with lines title "$\Phi(2u-2)$", "" using 1:3 with lines title "$\Phi(2u-1)$", "" using 1:4 with lines title "$\Phi(2u+0)$", "" using 1:5 with lines title "$\Phi(2u+1)$", "" using 1:6 with lines title "$\Phi(2u+2)$", "" using 1:7 with lines title "Sum of fine grids", "" using 1:8 with lines title "$\Phi(u)$"
#plot data_file using 1:2 with lines title '$\Phi(2u-2)$', "" using 1:3 with lines title '$\Phi(2u-1)$', "" using 1:4 with lines title '$\Phi(2u+0)$', "" using 1:5 with lines title '$\Phi(2u+1)$', "" using 1:6 with lines title '$\Phi(2u+2)$', "" using 1:7 with lines title "Sum of fine grids", "" using 1:8 with lines title '$\Phi(u)$'
#plot data_file using 1:2 with lines title '$\Phi(2u-2)$', "" using 1:3 with lines title '$\Phi(2u-1)$', "" using 1:4 with lines title '$\Phi(2u+0)$', "" using 1:5 with lines title '$\Phi(2u+1)$', "" using 1:6 with lines title '$\Phi(2u+2)$', "" using 1:7 with lines title '$\sum\limits_{v=-2}^{2}c_v\Phi(2u+v)$', "" using 1:8 with lines title '$\Phi(u)$'
#plot data_file using 1:2 with lines title '$\Phi(2u-2)$', "" using 1:3 with lines title '$\Phi(2u-1)$', "" using 1:4 with lines title '$\Phi(2u+0)$', "" using 1:5 with lines title '$\Phi(2u+1)$', "" using 1:6 with lines title '$\Phi(2u+2)$', "" using 1:7 with lines title '$\sum\limits_{v=-2}^{2}\Phi(2u+v)$', "" using 1:8 with lines title '$\Phi(u)$'
#plot data_file using 1:2 with lines title "{/Symbol F}(2u-2)", "" using 1:3 with lines title "{/Symbol F}(2u-1)", "" using 1:4 with lines title "{/Symbol F}(2u+0)", "" using 1:5 with lines title "{/Symbol F}(2u+1)", "" using 1:6 with lines title "{/Symbol F}(2u+2)", "" using 1:7 with lines title "Sum of fine grids", "" using 1:8 with lines title "{/Symbol F}(u)"
plot data_file using 1:2 with lines title "{/Symbol F}(2u-2)", "" using 1:3 with lines title "{/Symbol F}(2u-1)", "" using 1:4 with lines title "{/Symbol F}(2u+0)", "" using 1:5 with lines title "{/Symbol F}(2u+1)", "" using 1:6 with lines title "{/Symbol F}(2u+2)", "" using 1:7 with lines title "{/Symbol S}{/Symbol F}(2u+v)", "" using 1:8 with lines title "{/Symbol F}(u)"
set output
pause -1
quit
