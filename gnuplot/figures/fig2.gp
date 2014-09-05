reset
set term aqua enhanced
set termoption dashed
set xlabel 'u'
set ylabel 'f(u)'
set style data lines
set xrange [-5:5]
set yrange [-0.4:1.2]
plot data_file using 1:2 title "{/Symbol y}_3", "" using 1:3 title "{/Symbol y}_5", "" using 1:4 title "{/Symbol y}_7", "" using 1:5 title "{/Symbol y}_{/Symbol \245}"
pause -1
quit
