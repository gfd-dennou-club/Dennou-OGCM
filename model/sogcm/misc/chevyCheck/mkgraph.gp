dir="./../../unit_test/"

#
set terminal postscript eps enhanced color

#
set output '/dev/null'
plot for [i=1:8] dir.'waveFunc_'.i.'.dat' using 1:2 smooth cspline linewidth 4 title 'nMode='.i
set output "waveFunc.eps"
replot for [i=1:8] dir.'waveFunc_'.i.'.dat' using 1:2 w points title ''


#
set output '/dev/null'
plot for [i=1:4] dir.'expFunc_'.i.'.dat' using 1:2 smooth cspline linewidth 4 w l title 'lyrLenId='.i
replot for [i=5:6] dir.'expFunc_'.i.'.dat' using 1:2 linewidth 4 w l title 'lyrLenId='.i
set output "expFunc.eps"
replot for [i=1:6] dir.'expFunc_'.i.'.dat' using 1:2 w points title '' 


# comprison 
set output '/dev/null'
plot for [i=3:6] dir.'expFunc_'.i.'.dat' using 1:3 linewidth 4 w l title 'lyrLenId='.i
set output "expFunc_accuracy.eps"
replot for [i=3:6] dir.'expFunc_'.i.'.dat' using 1:4 linewidth 4 w l title 'lyrLenId='.i

# figure of error
set output "waveFunc_error.eps"
set logscale y
set format y "10^{%L}"
plot for [i=1:8] dir.'waveFunc_'.i.'.dat' using 1:5 w l linewidth 4 title 'nMode='.i
unset output

set output "expFunc_error.eps"
plot for [i=1:6] dir.'expFunc_'.i.'.dat' using 1:5 w l linewidth 4 title 'nLyrGPt='.i
unset output


