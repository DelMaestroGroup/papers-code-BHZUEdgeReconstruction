set terminal epslatex size 2.427,1.5 standalone color colortext 8
set border linewidth 3.25
set output 'Tot_N_vs_V0_for_12x3_GCE.tex'

set multiplot
#########################################################
set lmargin at screen 0.09
set rmargin at screen 0.985
set bmargin at screen 0.145
set tmargin at screen 0.98

set arrow nohead from 1.15,54 to 2.05,54 lw 2.5 lt 1 dt 3 lc rgb '#239b56'
set arrow nohead from 1.15,55 to 2.05,55 lw 2.5 lt 1 dt 3 lc rgb '#cb4335'
set arrow nohead from 1.15,56 to 2.05,56 lw 2.5 lt 1 dt 3 lc rgb '#2e86c1'

set ylabel '\fontsize{10}{10} $\langle N\rangle$' offset 6.4,-0.9
set xlabel '\fontsize{10}{10} $v_0$' offset 0,1.05
set xrange [1.15:1.75]
set yrange [53.5:56.2]

set ytics ('\fontsize{8}{8} $53$' 53, '\fontsize{8}{8} $54$' 54, '\fontsize{8}{8} $55$' 55, '\fontsize{8}{8} $56$' 56) offset 0.6,0

set xtics ('\fontsize{8}{8} $1.2$' 1.2,'\fontsize{8}{8} $1.4$' 1.4, '\fontsize{8}{8} $1.6$' 1.6, '\fontsize{8}{8} $1.7$' 1.7, '\fontsize{8}{8} $1.5$' 1.5,'\fontsize{8}{8} $1.3$' 1.3) offset 0,0.25

set label 1 '\scalebox{0.9}{$N_0$}' at 1.7,53.85 tc rgb '#239b56'
set label 2 '\scalebox{0.9}{$N_0+1$}' at 1.65,54.85 tc rgb '#cb4335'
set label 3 '\scalebox{0.9}{$N_0+2$}' at 1.65,55.85 tc rgb '#2e86c1'

plot "../../data/Tot_N_vs_V0_data_12x3_w4_GCE.txt" u (($1)/3):2 w lp lw 2.5 lt 1 pt 7 ps 1.2 lc rgb 'black' ti ''


#############################################################
unset multiplot
