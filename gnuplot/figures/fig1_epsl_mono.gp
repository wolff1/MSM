reset
#set term aqua
#set termoption dashed
set term epslatex mono
set output "bwfig1.tex"
set xlabel '$r$'
set ylabel '$g_l(r)$'
set label 2 '$(a)$' at graph 0.385, graph -0.10
set label 3 '$(2a)$' at graph 0.780, graph -0.10
set style data lines
set yrange [ 0.0 :   1.637 ]
set key nobox
#plot data_file using 1:2 title "g_0" ls 19, "" using 1:3 title "g_1" ls 15, "" using 1:4 title "g_2" ls 12, "" using 1:5 title "1/r" ls 16
#plot data_file using 1:2 title "g_0" ls 5, "" using 1:3 title "g_1" ls 3, "" using 1:4 title "g_2" ls 8, "" using 1:5 title "1/r" ls 4
#plot data_file using 1:2 title "g_0", "" using 1:3 title "g_1", "" using 1:4 title "g_2", "" using 1:5 title "1/r"
#plot data_file using 1:2 title "$g_0$", "" using 1:3 title "$g_1$", "" using 1:4 title "$g_2$", "" using 1:5 title "$\frac{1}{r}$"
#plot data_file using 1:2 title '$g_0$', "" using 1:3 title '$g_1$', "" using 1:4 title '$g_2$', "" using 1:5 title '$\frac{1}{r}$'
plot data_file using 1:2 title '$g_0$', "" using 1:3 title '$g_1$', "" using 1:4 title '$g_2$', "" using 1:5 title '$1/r$'
set output
pause -1
quit
