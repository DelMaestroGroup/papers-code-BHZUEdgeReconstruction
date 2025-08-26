set terminal epslatex size 2.427,1.5 standalone color colortext 8
set border linewidth 3.25
set output 'Tot_N_vs_V0_for_12x3_GCE.tex'

set multiplot
#########################################################
set lmargin at screen 0.11
set rmargin at screen 0.985
set bmargin at screen 0.125
set tmargin at screen 0.98

set arrow nohead from 1.15,54 to 1.75,54 lw 4 lt 1 dt 3 lc rgb '#239b56'
set arrow nohead from 1.15,55 to 1.75,55 lw 4 lt 1 dt 3 lc rgb '#cb4335'
set arrow nohead from 1.15,56 to 1.75,56 lw 4 lt 1 dt 3 lc rgb '#2e86c1'

set ylabel '\fontsize{10}{10} $\langle N\rangle$' offset 5.6,0.0
set xlabel '\fontsize{10}{10} $v_0$' offset 0,1.3
set xrange [1.15:1.75]
set yrange [53.5:56.2]

#set ytics ('\fontsize{8}{8} $53.0$' 53.0, '\fontsize{8}{8} $53.5$' 53.5, '\fontsize{8}{8} $54.0$' 54.0, '\fontsize{8}{8} $54.5$' 54.5, '\fontsize{8}{8} $55.0$' 55.0, '\fontsize{8}{8} $55.5$' 55.5, '\fontsize{8}{8} $56.0$' 56) offset 0.6,0

set ytics ('\fontsize{8}{8} $53$' 53, '\fontsize{8}{8} $54$' 54, '\fontsize{8}{8} $55$' 55, '\fontsize{8}{8} $56$' 56) offset 0.9,0

#set xtics ('\fontsize{8}{8} $3.5$' 3.5,'\fontsize{8}{8} $4.0$' 4.0, '\fontsize{8}{8} $4.5$' 4.5, '\fontsize{8}{8} $5.0$' 5.0, '\fontsize{8}{8} $5.5$' 5.5, '\fontsize{8}{8} $6.0$' 6.0,'\fontsize{8}{8} $6.5$' 6.5) offset 0,0.25

#set xtics ('\fontsize{8}{8} $1.2$' 1.2,'\fontsize{8}{8} $1.4$' 1.4, '\fontsize{8}{8} $1.6$' 1.6, '\fontsize{8}{8} $1.8$' 1.8, '\fontsize{8}{8} $2.0$' 2.0,'\fontsize{8}{8} $6.5$' 6.5) offset 0,0.25

set xtics ('\fontsize{8}{8} $1.2$' 1.2,'\fontsize{8}{8} $1.4$' 1.4, '\fontsize{8}{8} $1.6$' 1.6, '\fontsize{8}{8} $1.7$' 1.7, '\fontsize{8}{8} $1.5$' 1.5,'\fontsize{8}{8} $1.3$' 1.3) offset -0.3,0.5

#set label 1 '\scalebox{0.9}{$N_0$}' at 1.95,54.15 tc rgb '#239b56'
#set label 2 '\scalebox{0.9}{$N_0+1$}' at 1.9,55.15 tc rgb '#cb4335'
#set label 3 '\scalebox{0.9}{$N_0+2$}' at 1.9,55.85 tc rgb '#2e86c1'

set label 1 '\scalebox{0.9}{$N_0$}' at 1.7,53.85 tc rgb '#239b56'
set label 2 '\scalebox{0.9}{$N_0+1$}' at 1.65,54.85 tc rgb '#cb4335'
set label 3 '\scalebox{0.9}{$N_0+2$}' at 1.65,55.85 tc rgb '#2e86c1'

plot "../../data/Tot_N_vs_V0_data_12x3_w4_GCE.txt" u (($1)/3):2 w lp lw 4 lt 1 pt 7 ps 1.2 lc rgb 'black' ti ''

#plot '-' u 1:2 w p pt 7 ps 1.6 lc rgb '#239b56' ti '','-' u 1:2 w p pt 7 ps 1.6 lc rgb '#cb4335' ti '','-' u 1:2 w p pt 7 ps 1.6 lc rgb '#2e86c1' ti '', "Tot_N_vs_V0_data_12x3_w4_GCE.txt" u 1:2 w lp lw 4 lt 1 pt 7 ps 1.2 lc rgb 'black' ti ''

#N0 point: (4.6, 54)
#4.6 54
#e

#N0+1 point: (4.2, 55)
#4.2 55.206
#e

#N0+2 point: (4.0, 56)
#4.0 56
#e


#############################################################
unset multiplot
