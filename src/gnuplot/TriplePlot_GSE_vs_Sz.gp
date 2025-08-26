set terminal epslatex size 2.1,2 standalone color colortext 8
set output 'GSE_vs_Sz_for_diff_Vo_12x3_triple.tex'

set border behind linewidth 2.75

set multiplot
#########################################################
set lmargin at screen 0.15
set rmargin at screen 0.99
set tmargin at screen 0.98
set bmargin at screen 0.71

set xtics nomirror

set yr [-0.002:0.022]
set xr [-2.2:2.2]
set xtics ('' -2, '' -1.5, '' -1, '' -0.5, '' 0, '' 0.5, '' 1, '' 1.5, '' 2) offset -0.2,0.25

set ytics rotate by 0 offset 0.7,0.1 ('\fontsize{7}{7} $0.0$' 0.000,'\fontsize{7}{7} $0.5$' 0.005, '\fontsize{7}{7} $1.0$' 0.010, '\fontsize{7}{7} $1.5$' 0.015, '\fontsize{7}{7} $2.0$' 0.020)

set key width -3 samplen 0.1 spacing 0.7 at 0.9,0.0211
set label 1 '\scalebox{1.0}{(a)}' at -2.1,0.019
set label 2 '\scalebox{0.85}{$N_0$}' at 1.6,0.0195

plot "../../data/GSE_vs_Sz_for_diff_Vo4.6.txt" u 1:((($2)+88.7462080)/1.0) w l lw 2.5 lt 1 dt 3 lc rgb '#abebc6' ti '',"" u 1:((($2)+88.7462080)/1.0) w p ps 1.0 pt 7 lc rgb '#239b56' ti '\scalebox{0.9}{$v_0=1.53$}', "" using ($1 -0.07):((($2)+88.7462080)/1.0):(0.14):(0) w vectors nohead lt 1 lw 3 lc rgb '#239b56' ti ''
unset ytics
unset key
unset label 1
unset label 2
#########################################################

set lmargin at screen 0.15
set rmargin at screen 0.99
set tmargin at screen 0.68
set bmargin at screen 0.41

set yr [-0.0015:0.015]
set ylabel '\fontsize{8}{8} $[E(S^z_{\mathrm{tot}})-E_{GS}]\times 10^{-2}$' offset 5.95,0

set ytics rotate by 0 offset 0.7,0.15 ('\fontsize{7}{7} $0.0$' 0.000, '\fontsize{7}{7} $0.4$' 0.004, '\fontsize{7}{7} $0.8$' 0.008, '\fontsize{7}{7} $1.2$' 0.012)

set key width -2 samplen 0.1 spacing 0.7 at 0.9,0.0144
set label 1 '\scalebox{1.0}{(b)}' at -2.1,0.013
set label 2 '\scalebox{0.85}{$N_0+1$}' at 1.4,0.0065

plot "../../data/GSE_vs_Sz_for_diff_Vo4.2.txt" u 1:((($2)+88.8977663)/1.0) w l lw 2.5 lt 1 dt 3 lc rgb  '#f1948a' ti '', "" u 1:((($2)+88.8977663)/1.0) w p ps 1.0 pt 7 lc rgb '#cb4335' ti '\scalebox{0.9}{$v_0=1.40$}',"" using ($1 -0.07):((($2)+88.8977663)/1.0):(0.14):(0) w vectors nohead lt 1 lw 3 lc rgb '#cb4335' ti ''

unset ylabel
unset key
unset ytics
unset label 1
unset label 2
#########################################################

set lmargin at screen 0.15
set rmargin at screen 0.99
set tmargin at screen 0.38
set bmargin at screen 0.11

set yr [-0.001:0.010]
set xlabel '\fontsize{8}{8} $S^{z}_{\mathrm{tot}}$' offset 0.0,1.1

set xtics ('\fontsize{7}{7} $-2$' -2, '\fontsize{7}{7} $-\frac{3}{2}$' -1.5, '\fontsize{7}{7} $-1$' -1, '\fontsize{7}{7} $-\frac{1}{2}$' -0.5, '\fontsize{7}{7} $0$' 0, '\fontsize{7}{7} $\frac{1}{2}$' 0.5, '\fontsize{7}{7} $1$' 1, '\fontsize{7}{7} $\frac{3}{2}$' 1.5, '\fontsize{7}{7} $2$' 2) offset -0.2,0.4

set ytics rotate by 0 offset 0.7,-0.05 ('\fontsize{7}{7} $0.0$' 0.000,'\fontsize{7}{7} $0.3$' 0.003, '\fontsize{7}{7} $0.6$' 0.006, '\fontsize{7}{7} $0.9$' 0.009)

set key width -2 samplen 0.1 spacing 0.7 at 0.85,0.0019
set label 1 '\scalebox{1.0}{(c)}' at -2.1,0.0013
set label 2 '\scalebox{0.85}{$N_0+2$}' at 1.4,0.0013

plot "../../data/GSE_vs_Sz_for_diff_Vo4.0.txt" u 1:((($2)+89.0718564)/1.0) w l lw 2.5 lt 1 dt 3 lc rgb '#aed6f1' ti '', "" u 1:((($2)+89.0718564)/1.0) w p ps 1.0 pt 7 lc rgb '#2e86c1' ti '\scalebox{0.9}{$v_0=1.33$}',"" using ($1 -0.07):((($2)+89.0718564)/1.0):(0.14):(0) w vectors nohead lt 1 lw 3 lc rgb '#2e86c1' ti ''

#############################################################
unset multiplot
