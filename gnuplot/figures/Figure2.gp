# Figure 2 - Sinc test
reset
#set term windows mono
#set term epslatex mono
set term postscript eps enhanced mono font ',21'
set output "bwfig2.eps"
set xlabel 'u'
set ylabel 'f(u)'
set style data lines
set key nobox
set xrange [-5:5]
set yrange [-0.4:1.2]
plot "../../data/figures/Figure2.dat" every :::3::3 using ($2):($3) title '{/Symbol y}_{/Symbol \245}', \
"" every :::0::0 using ($2):($3) title '{/Symbol y}_3', \
"" every :::1::1 using ($2):($3) title '{/Symbol y}_5', \
"" every :::2::2 using ($2):($3) title '{/Symbol y}_7'
set output
#pause -1
quit