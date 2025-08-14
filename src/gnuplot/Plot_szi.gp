set terminal epslatex size 2.1,2 standalone color colortext 8
set output 'szi_12x3_w4.tex'

set border behind linewidth 2.75

set multiplot
#########################################################
set xr [-0.5:11.5]
set xtics ('' 0, '' 2, '' 4, '' 6, '' 8, '' 10, '' 12) offset -0.3,0.25

set yr [-0.07:0.30]
set ytics ('\fontsize{7}{7} $-0.1$' -0.1, '\fontsize{7}{7} $0.0$' 0.0, '\fontsize{7}{7} $0.1$' 0.1, '\fontsize{7}{7} $0.2$' 0.2, '\fontsize{7}{7} $0.3$' 0.3) offset 0.7,-0.1

set xzeroaxis lw 2.5 lt 1 dt 3 lc rgb 'black'
set arrow nohead from 8,-0.07 to 8,0.3 lw 2.5 lt 1 dt 3 lc rgb 'black'

set lmargin at screen 0.14
set rmargin at screen 0.972
set tmargin at screen 0.98
set bmargin at screen 0.71

set key width -3 samplen 0.1 spacing 0.7 at 11.5,0.28
set label 1 '\scalebox{1.0}{(a)} \scalebox{0.8}{$N_0$, $S^z_{tot}=0$}' at 0.1,0.24

plot "../data/Sz_rx_yavg_Nup27_Ndn27_L12_W03_V4.60_Sz0.0.txt" u (($1)-1.0):2 w lp ps 0.9 pt 7 lw 4.5 lt 1 lc rgb '#cb4335' ti '\fontsize{7}{7} $\alpha=s$', "" u (($1)-1.0):3 w lp ps 0.9 pt 5 lw 4.5 lt 1 lc rgb '#2e86c1' ti '\fontsize{7}{7} $\alpha=p$'

unset label 1
unset key
#########################################################

set lmargin at screen 0.14
set rmargin at screen 0.972
set tmargin at screen 0.68
set bmargin at screen 0.41

set ylabel '\fontsize{8}{8} $\langle S^z_{r_x,\alpha}\rangle$' offset 6.1,0
set label 1 '\scalebox{1.0}{(b)} \scalebox{0.8}{$N_0+1$, $S^z_{tot}=1/2$}' at 0.1,0.24

plot "../data/Sz_rx_yavg_Nup28_Ndn27_L12_W03_V4.20_Sz0.5.txt" u (($1)-1.0):2 w lp ps 0.9 pt 7 lw 4.5 lt 1 lc rgb '#cb4335' ti '', "" u (($1)-1.0):3 w lp ps 0.9 pt 5 lw 4.5 lt 1 lc rgb '#2e86c1' ti ''

unset ylabel
unset label 1
#########################################################

set lmargin at screen 0.14
set rmargin at screen 0.972
set tmargin at screen 0.38
set bmargin at screen 0.11

set xlabel '\fontsize{8}{8} $r_x$' offset 0.0,1.15
set xtics ('\fontsize{7}{7} $0$' 0, '\fontsize{7}{7} $2$' 2, '\fontsize{7}{7} $4$' 4, '\fontsize{7}{7} $6$' 6, '\fontsize{7}{7} $8$' 8, '\fontsize{7}{7} $10$' 10, '\fontsize{7}{7} $12$' 12) offset -0.3,0.4

set label 1 '\scalebox{1.0}{(c)} \scalebox{0.8}{$N_0+2$, $S^z_{tot}=1$}' at 0.1,0.24

plot "../data/Sz_rx_yavg_Nup29_Ndn27_L12_W03_V4.00_Sz1.0.txt" u (($1)-1.0):2 w lp ps 0.9 pt 7 lw 4.5 lt 1 lc rgb '#cb4335' ti '', "" u (($1)-1.0):3 w lp ps 0.9 pt 5 lw 4.5 lt 1 lc rgb '#2e86c1' ti ''

#############################################################
unset multiplot
