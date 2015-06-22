# Figure 1 - Splitting test
set term windows
#set terminal postscript eps color
#set size 0.75,0.75
set xlabel "r"
set ylabel "g_l(r)"
set label 2 "(a)" at graph 0.39, graph -0.08
set label 3 "(2a)" at graph 0.785, graph -0.08
set title 'Kernel Splitting for 3-level MSM'
#set grid
set style data linespoints
set yrange [ 0.0 : 2.0 ]
set key box
plot "Figure1.dat" every :::0::0 using ($2):($3) title "g_0" lc rgb "black", \
"Figure1.dat" every :::1::1 using ($2):($3) title "g_1" lc rgb "black", \
"Figure1.dat" every :::2::2 using ($2):($3) title "g_2" lc rgb "black", \
"Figure1.dat" every :::3::3 using ($2):($3) title "1/r" lc rgb "black"


#using 1:%d title \"g_%d\" lc rgb \"black\",
#data_file using 1:%d title \"1/r\" lc rgb \"black\"

pause -1
quit