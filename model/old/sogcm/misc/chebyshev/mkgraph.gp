dir="./data/"
k0="20"
k1="40"
k2="60"
k3="80"

lyr1="0.05"
lyr2="0.01"
lyr3="0.005"
lyr4="0.0025"
lyr5="0.001"

#
set terminal postscript eps enhanced color


do for [k=0:3] {
   kMax = value(sprintf('k%i',k))

   #
   set output '/dev/null'
   set size 1.35, 1.0
   set yrange [-1:1]
   plot for [i=1:8] dir.'waveFunc_'.i.'_k'.kMax.'.dat' using 1:2 smooth cspline linewidth 4 title 'mode='.i
   set output dir.'waveFunc_dist_k'.kMax.'.eps'
   replot for [i=1:8] dir.'waveFunc_'.i.'_k'.kMax.'.dat' using 1:7 w points pt 19 title ''


   #
   set output '/dev/null'
   set yrange [-0.05:1]
   plot for [i=1:2] dir.'expFunc_'.i.'_k'.kMax.'.dat' using 1:2 smooth cspline linewidth 4 w l title 'lyrLen='.value(sprintf("lyr%i",i))
   replot for [i=3:5] dir.'expFunc_'.i.'_k'.kMax.'.dat' using 1:2 linewidth 4 w l title 'lyrLen='.value(sprintf("lyr%i",i))
   set output dir."expFunc_dist_k".kMax.".eps"
   replot for [i=1:5] dir.'expFunc_'.i.'_k'.kMax.'.dat' using 1:7 w points pt 19  title '' 


   # comprison 
   set output '/dev/null'
   set yrange [-0.01:1.2]
   plot for [i=1:5] dir.'expFunc_'.i.'_k'.kMax.'.dat' using 1:3 linewidth 4 w l title 'lyrLen='.value(sprintf("lyr%i",i))
   set output dir."expFunc_accuracy_k".kMax.".eps"
   replot for [i=1:5] dir.'expFunc_'.i.'_k'.kMax.'.dat' using 1:4 linewidth 5 lt 0 w l title 'answer lyrLen='.value(sprintf("lyr%i",i))

   # figure of error
   set output '/dev/null'
   set yrange [1e-18:1e0]
   set output dir."waveFunc_error_k".kMax.".eps"
   set logscale y
   set format y "10^{%L}"
   plot for [i=1:8] dir.'waveFunc_'.i.'_k'.kMax.'.dat' using 1:5 w l linewidth 4 title 'mode='.i
   unset output

   set output dir."expFunc_error_k".kMax.".eps"
   plot for [i=1:5] dir.'expFunc_'.i.'_k'.kMax.'.dat' using 1:5 w l linewidth 4 title 'lyrLen='.value(sprintf("lyr%i",i))
   unset output

   unset logscale
   unset format
}
