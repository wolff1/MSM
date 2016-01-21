# 1D interpolation data
reset
set term windows
set xlabel 's'
set ylabel 'f(s)'
set style data lines
set key box

# USE MULTI-PLOT
set multiplot layout 3,1 rowsfirst

# First plot - gamma vs interpolated gamma
plot "1DInterp.dat" using ($2):($3) title "Gamma(s)", \
"" using ($2):($4) title "AGamma(s)"

# Second plot - abs error
plot "" using ($2):($5) title "Abs Diff"

# Third plot - rel error
plot "" using ($2):($6) title "Rel Diff"

unset multiplot

set output
pause -1
quit