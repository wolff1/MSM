# Figure 1 - Splitting test
reset
#set term windows mono
#set term epslatex mono
set term postscript eps enhanced mono
set output "bwfig1.eps"
set xlabel 'r'
set ylabel 'g_l(r)'
set label 2 '(a)' at graph 0.385, graph -0.10
set label 3 '(2a)' at graph 0.780, graph -0.10
set style data lines
set yrange [ 0.0 : 1.637 ]
set key box
plot "../../data/figures/Figure1.dat" every :::0::0 using ($2):($3) title "g_0", \
"" every :::1::1 using ($2):($3) title "g_1", \
"" every :::2::2 using ($2):($3) title "g_2", \
"" every :::3::3 using ($2):($3) title "1/r"
set output
#pause -1
quit