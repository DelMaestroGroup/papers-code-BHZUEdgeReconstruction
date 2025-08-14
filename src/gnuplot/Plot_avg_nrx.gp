set terminal epslatex size 4.0,1.5 standalone color colortext 8
set border linewidth 3.25
set output 'Avg_nrx_vs_rx_for_12x3_DMRG.tex'

set multiplot
#########################################################
set ylabel '\fontsize{10}{10} $\langle n_{r_x,\alpha}\rangle$' offset 5.3,0.0
set xlabel '\fontsize{10}{10} $r_x$' offset 0,1.1

set xrange [-0.5:11.5]
set yrange [0.0:2.0]

set ytics nomirror
set xtics nomirror

set ytics ('\fontsize{8}{8} $0.0$' 0.0, '\fontsize{8}{8} $0.4$' 0.4, '\fontsize{8}{8} $0.8$' 0.8, '\fontsize{8}{8} $1.2$' 1.2, '\fontsize{8}{8} $1.6$' 1.6, '\fontsize{8}{8} $2.0$' 2.0) offset 0.6,-0.1

set xtics ('\fontsize{7}{7} $0$' 0, '\fontsize{7}{7} $2$' 2, '\fontsize{7}{7} $4$' 4, '\fontsize{7}{7} $6$' 6, '\fontsize{7}{7} $8$' 8, '\fontsize{7}{7} $10$' 10, '\fontsize{7}{7} $12$' 12) offset 0.0,0.4

set arrow nohead from 8,0.0 to 8,2.0 lw 2.5 lt 1 dt 3 lc rgb 'black'

set lmargin at screen 0.09
set rmargin at screen 0.38
set bmargin at screen 0.14
set tmargin at screen 0.975

set key width -3 samplen 0.1 spacing 0.7 at 11.5,1.95
set label 1 '\scalebox{1.0}{(a)} \scalebox{0.9}{$v_0=1.53$}' at -0.2,1.85

plot "../../data/raw_n_rx_L12_W03_V4.60_Nup27_Ndn27.txt" u (($1)-1.0):(($2)/3.0) w lp ps 0.9 pt 7 lw 2.5 lt 1 lc rgb '#cb4335' ti '\fontsize{7}{7} $s$', "" u (($1)-1.0):(($3)/3.0) w lp ps 0.9 pt 5 lw 2.5 lt 1 lc rgb '#2e86c1' ti '\fontsize{7}{7} $p$'

unset ylabel
unset label 1
unset key

set lmargin at screen 0.39
set rmargin at screen 0.68
set bmargin at screen 0.14
set tmargin at screen 0.975

set ytics ('' 0.0, '' 0.4, '' 0.8, '' 1.2, '' 1.6, '' 2.0) offset 0.6,0
set label 1 '\scalebox{1.0}{(b)} \scalebox{0.9}{$v_0=1.40$}' at -0.2,1.85

plot "../../data/raw_n_rx_L12_W03_V4.20_Nup28_Ndn27.txt" u (($1)-1.0):(($2)/3.0) w lp ps 0.9 pt 7 lw 2.5 lt 1 lc rgb '#cb4335' ti '', "" u (($1)-1.0):(($3)/3.0) w lp ps 0.9 pt 5 lw 2.5 lt 1 lc rgb '#2e86c1' ti ''

unset label 1

set lmargin at screen 0.69
set rmargin at screen 0.98
set bmargin at screen 0.14
set tmargin at screen 0.975

set label 1 '\scalebox{1.0}{(c)} \scalebox{0.9}{$v_0=1.33$}' at -0.2,1.85

plot "../../data/raw_n_rx_L12_W03_V4.00_Nup29_Ndn27.txt" u (($1)-1.0):(($2)/3.0) w lp ps 0.9 pt 7 lw 2.5 lt 1 lc rgb '#cb4335' ti '', "" u (($1)-1.0):(($3)/3.0) w lp ps 0.9 pt 5 lw 2.5 lt 1 lc rgb '#2e86c1' ti ''

#############################################################
unset multiplot
