reset
#set term aqua
#set termoption dashed
set terminal epslatex standalone color
set output "fig2.tex"
set xlabel 'u'
set ylabel 'f(u)'
set style data lines
set key box
set xrange [-5:5]
set yrange [-0.4:1.2]
#plot data_file using 1:2 title "Cubic" ls 12, "" using 1:3 title "Quintic" ls 19, "" using 1:4 title "Septic" ls 15, "" using 1:5 title "Sinc" ls -1
#plot data_file using 1:2 title "Cubic" ls 4, "" using 1:3 title "Quintic" ls 7, "" using 1:4 title "Septic" ls 8, "" using 1:5 title "Sinc" ls -1
plot data_file using 1:2 title "Cubic", "" using 1:3 title "Quintic", "" using 1:4 title "Septic", "" using 1:5 title "Sinc"
set output
pause -1
quit
