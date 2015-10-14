# dgamma1.gp - visualize gamma derivatives over an interval
reset
set term windows mono
#set term epslatex mono
#set term postscript eps enhanced mono font ',21'
#set output "dgamma1.eps"
set title 'i-th derivative of gamma'
set xlabel 's'
set ylabel 'g^{(i)}'
set style data lines
set key box
#plot "../data/dgamma_2015-10-14_11-40-21.dat" every :::0::0 using ($2):($3) title "i=0", \
#										   "" every :::1::1 using ($2):($3) title "i=1", \
#										   "" every :::2::2 using ($2):($3) title "i=2", \
#										   "" every :::3::3 using ($2):($3) title "i=3", \
#										   "" every :::4::4 using ($2):($3) title "i=4"
plot "../data/dgamma_2015-10-14_11-53-57.dat" every :::4::4 using ($2):($3) title "i=4"
set output
pause -1
quit