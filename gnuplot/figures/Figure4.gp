# Figure 4 - B-spline nesting
reset
#set term windows mono
#set term epslatex mono
set term postscript eps enhanced mono
set output "bwfig4.eps"
set xlabel 'u'
set ylabel 'f(u)'
set style data lines
set yrange [   0.000 :   0.833 ]
#set key at 0.10,0.50 bottom left nobox
plot "../../data/figures/Figure4.dat" using ($2):($3) title '{/Symbol F}(2u-2)', \
"" using ($2):($4) title '{/Symbol F}(2u-1)', \
"" using ($2):($5) title '{/Symbol F}(2u+0)', \
"" using ($2):($6) title '{/Symbol F}(2u+1)', \
"" using ($2):($7) title '{/Symbol F}(2u+2)', \
"" using ($2):($8) title '{/Symbol S}{/Symbol F}(2u+v)', \
"" using ($2):($9) title '{/Symbol F}(u)'
set output
#pause -1
quit