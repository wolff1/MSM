# tmp_data_template.gp (Created 07/31/2014 - MAW)
#
# Usage:
#  gnuplot < tmp_data_template.gp
#
#set term pngcairo
#set output "damped_sine.png"
#set term windows
set term aqua
set xlabel 'x'
set ylabel 'f(x)'
#set title 'The title of the plot'
set grid
set style data lines
plot data_file using 1:2 lw 3 linecolor rgb "blue"
pause -1
quit
