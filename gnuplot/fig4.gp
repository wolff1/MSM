set term aqua
set termoption dash
set xlabel 'u'
set ylabel 'f(u)'
set title 'Nesting of Spline Function Spaces'
set grid
set style data lines
set yrange [   0.000 :   0.833 ]
set key box
plot data_file using 1:2 with lines title "Phi(2u-2)",data_file using 1:3 with lines title "Phi(2u-1)",data_file using 1:4 with lines title "Phi(2u+0)",data_file using 1:5 with lines title "Phi(2u+1)",data_file using 1:6 with lines title "Phi(2u+2)",data_file using 1:7 with lines title "Sum of fine grids",data_file using 1:8 with lines title "Phi(u)"
pause -1
quit
