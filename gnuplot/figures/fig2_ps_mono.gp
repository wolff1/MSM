reset
#set term aqua
#set termoption dashed
#set term epslatex mono
set term postscript eps enhanced mono
set output "bwfig2.eps"
#set xlabel '$u$'
#set ylabel '$f(u)$'
set xlabel 'u'
set ylabel 'f(u)'
set style data lines
set key nobox
set xrange [-5:5]
set yrange [-0.4:1.2]
#plot data_file using 1:2 title "Cubic" ls 12, "" using 1:3 title "Quintic" ls 19, "" using 1:4 title "Septic" ls 15, "" using 1:5 title "Sinc" ls -1
#plot data_file using 1:2 title "Cubic" ls 4, "" using 1:3 title "Quintic" ls 7, "" using 1:4 title "Septic" ls 8, "" using 1:5 title "Sinc" ls -1
#plot data_file using 1:2 title "Cubic", "" using 1:3 title "Quintic", "" using 1:4 title "Septic", "" using 1:5 title "Sinc"
#plot data_file using 1:2 title "$\Phi_3$", "" using 1:3 title "$\Phi_5$", "" using 1:4 title "$\Phi_7$", "" using 1:5 title "$\Psi_\infnty$"
#plot data_file using 1:2 title '$\Psi_3$', "" using 1:3 title '$\Psi_5$', "" using 1:4 title '$\Psi_7$', "" using 1:5 title '$\Psi_{\infty}$'
plot data_file using 1:2 title "{/Symbol y}_3", "" using 1:3 title "{/Symbol y}_5", "" using 1:4 title "{/Symbol y}_7", "" using 1:5 title "{/Symbol y}_{/Symbol \245}"
set output
pause -1
quit
