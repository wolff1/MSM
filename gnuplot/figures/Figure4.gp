# Figure 4 - B-spline nesting
set term windows
#set term windows mono
#set terminal postscript eps color
set xlabel "u"
set ylabel "f(u)"
#set title 'Nesting of Spline Function Spaces'
#set grid
set style data lines
set key box
plot "Figure4.dat" using ($2):($3) title 'Phi-2', \
"Figure4.dat" using ($2):($4) title 'Phi-1', \
"Figure4.dat" using ($2):($5) title 'Phi-0', \
"Figure4.dat" using ($2):($6) title 'Phi+1', \
"Figure4.dat" using ($2):($7) title 'Phi+2', \
"Figure4.dat" using ($2):($8) title 'Sum(Phi+i)', \
"Figure4.dat" using ($2):($9) title 'CoarsePhi'

pause -1
quit