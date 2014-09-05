reset
set term aqua enhanced
set termoption dashed
set xlabel 'r'
set ylabel 'g_l(r)'
set label 2 '(a)' at graph 0.385, graph -0.07
set label 3 '(2a)' at graph 0.780, graph -0.07
set style data lines
set yrange [ 0.0 :   1.637 ]
plot data_file using 1:2 title "g_0", "" using 1:3 title "g_1", "" using 1:4 title "g_2", "" using 1:5 title "1/r"
pause -1
quit
