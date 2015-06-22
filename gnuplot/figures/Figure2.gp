# Figure 1 - Splitting test
set term windows mono
#set terminal postscript eps color
set xlabel "u"
set ylabel "f(u)"
set title 'Fundamental splines of varying degree and sinc function'
#set grid
set style data lines
set key box
set xrange[-5:5]
plot "Figure2.dat" every :::0::0 using ($2):($3) title '$\Psi_3$', \
"Figure2.dat" every :::1::1 using ($2):($3) title '$\Psi_5$', \
"Figure2.dat" every :::2::2 using ($2):($3) title '$\Psi_7$', \
"Figure2.dat" every :::3::3 using ($2):($3) title '$\Psi_{\infty}$'

pause -1
quit