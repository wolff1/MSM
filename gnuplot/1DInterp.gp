# 1D interpolation data
reset
set term windows
#set term epslatex mono
#set term postscript eps enhanced mono font ',21'
#set term postscript eps enhanced font ',21'
#set output "gamma_deriv.eps"
set xlabel 's'
set ylabel 'f(s)'
set style data lines
set key box
#set size ratio -1
#set yrange [0 : 2.5]

#plot "1DInterp.dat" using ($2):($3) title "Gamma(s)", \
#"" using ($2):($4) title "AGamma(s)", \
#"" using ($2):($5) title "Abs Diff", \
#"" using ($2):($6) title "Rel Diff"

# USE MULTI-PLOT
set multiplot layout 2,1 rowsfirst

# First plot - gamma vs interpolated gamma
plot "1DInterp.dat" using ($2):($3) title "Gamma(s)", \
"" using ($2):($4) title "AGamma(s)"

# Second plot - abs and rel difference
plot "1DInterp.dat" using ($2):($5) title "Abs Diff", \
"" using ($2):($6) title "Rel Diff"

unset multiplot

set output
pause -1
quit