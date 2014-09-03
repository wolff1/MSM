set term aqua
set xlabel 'r'
set ylabel 'g_l(r)'
set label 2 '(a)' at graph 0.39, graph -0.08
set label 3 '(2a)' at graph 0.785, graph -0.08
set title 'Kernel Splitting for 3-level MSM'
set grid
set style data linespoints
set yrange [ 0.0 :   1.637 ]
set key box
plot data_file using 1:2 title "g_0" lc rgb "black",data_file using 1:3 title "g_1" lc rgb "black",data_file using 1:4 title "g_2" lc rgb "black",data_file using 1:5 title "1/r" lc rgb "black"
pause -1
quit
