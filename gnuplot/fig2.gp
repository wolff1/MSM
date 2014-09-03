set term aqua
set xlabel 'u'
set ylabel 'f(u)'
set title 'Fundamental splines of varying degree and sinc function'
set style data lines
set key box
plot [][-0.4:1.2] 0.0 title "", data_file using 1:2 title "Cubic", "" using 1:3 title "Quintic", "" using 1:4 title "Septic", "" using 1:5 title "Sinc"
pause -1
quit
