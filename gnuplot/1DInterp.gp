# 1D interpolation data
reset
set term windows
set xlabel 's'
set ylabel 'f(s)'
set style data lines
set key box

stats "1DInterp.dat" using 2:5 prefix "A" nooutput
stats "" using 2:6 prefix "R" nooutput

# USE MULTI-PLOT
set multiplot layout 3,1 rowsfirst

# First plot - gamma vs interpolated gamma
set object circle at 1.0,1.0 size 0.05
set label 1  at R_pos_max_y,0.0 "+"
plot "" using ($2):($3) title "Gamma(s)", \
"" using ($2):($4) title "AGamma(s)"
unset object
unset label

# Second plot - abs error
set object circle at 1.0,(A_max_y-A_min_y)/2.0 size 0.05
set label 2 at R_pos_max_y,0.0 "+"
plot "" using ($2):($5) title "Abs Diff"
unset object
unset label

# Third plot - rel error
set object circle at 1.0,(R_max_y-R_min_y)/2.0 size 0.05
set label 3 at R_pos_max_y,0.0 "+"
plot "" using ($2):($6) title "Rel Diff"
unset object
unset label

unset multiplot

set output
pause -1
quit