set terminal epslatex size 2.1,2 standalone color colortext 8
set output 'del_ni_12x3_w4.tex'

set border behind linewidth 2.75

set multiplot
#########################################################
set xr [-0.5:11.5]
set yr [-0.01:0.55]

set ytics ('\fontsize{7}{7} $0.0$' 0.0, '\fontsize{7}{7} $0.1$' 0.1, '\fontsize{7}{7} $0.2$' 0.2, '\fontsize{7}{7} $0.3$' 0.3, '\fontsize{7}{7} $0.4$' 0.4, '\fontsize{7}{7} $0.5$' 0.5) offset 0.7,-0.1

set lmargin at screen 0.137
set rmargin at screen 0.972
set tmargin at screen 0.985
set bmargin at screen 0.555

set ylabel '\fontsize{8}{8} $\Delta n$' offset 6.1,0
set xtics ('' 0, '' 2, '' 4, '' 6, '' 8, '' 10, '' 12) offset -0.3,0.25

set key width -3 samplen 0.1 spacing 0.7 at 11.5,0.5
set label 1 '\scalebox{1.0}{(a)} \scalebox{0.8}{$\langle n_{r_x,\alpha}^{N_0+1}\rangle - \langle n_{r_x,\alpha}^{N_0}\rangle$}' at 0.1,0.45

set arrow nohead from 8,-0.01 to 8,0.55 lw 2.5 lt 1 dt 3 lc rgb 'black'

plot "../data/delta_n_rx_L12_W03_V1_4.20_V2_4.60.txt" u (($1)-1.0):(($2)/3.0) w lp ps 0.9 pt 7 lw 4.5 lt 1 lc rgb '#cb4335' ti '\fontsize{7}{7} $\alpha=s$', "" u (($1)-1.0):(($3)/3.0) w lp ps 0.9 pt 5 lw 4.5 lt 1 lc rgb '#2e86c1' ti '\fontsize{7}{7} $\alpha=p$'

unset key
unset label 1
unset key
#########################################################

set lmargin at screen 0.137
set rmargin at screen 0.972
set tmargin at screen 0.525
set bmargin at screen 0.095

set xlabel '\fontsize{8}{8} $r_x$' offset 0.0,1.2
set xtics ('\fontsize{7}{7} $0$' 0, '\fontsize{7}{7} $2$' 2, '\fontsize{7}{7} $4$' 4, '\fontsize{7}{7} $6$' 6, '\fontsize{7}{7} $8$' 8, '\fontsize{7}{7} $10$' 10, '\fontsize{7}{7} $12$' 12) offset -0.3,0.4

set label 1 '\scalebox{1.0}{(b)} \scalebox{0.8}{$\langle n_{r_x,\alpha}^{N_0+2}\rangle - \langle n_{r_x,\alpha}^{N_0}\rangle$}' at 0.1,0.45

plot "../data/delta_n_rx_L12_W03_V1_4.00_V2_4.60.txt" u (($1)-1.0):(($2)/3.0) w lp ps 0.9 pt 7 lw 4.5 lt 1 lc rgb '#cb4335' ti '', "" u (($1)-1.0):(($3)/3.0) w lp ps 0.9 pt 5 lw 4.5 lt 1 lc rgb '#2e86c1' ti ''

unset ylabel
unset key
unset ytics
unset label 1
#############################################################
unset multiplot
