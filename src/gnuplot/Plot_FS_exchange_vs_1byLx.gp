set terminal epslatex size 2.427,1.5 standalone color colortext 8
set border linewidth 3.25
set output 'FS_exchange_vs_1byLx_Nx3_DMRG.tex'

set multiplot
#########################################################
set lmargin at screen 0.13
set rmargin at screen 0.965
set bmargin at screen 0.155
set tmargin at screen 0.99

set ylabel '\fontsize{10}{10} $J\times 10^{-2}$' offset 6.1,0
set xlabel '\fontsize{10}{10} $1/L_x$' offset 0.4,1.0
set xrange [0.0:0.1]
set yrange [*:*]

f(x) = 0.034689635 - 0.814004893*(x - 0.0625)

#set ytics ('\fontsize{8}{8} $1.6$' 0.016, '\fontsize{8}{8} $2.0$' 0.020, '\fontsize{8}{8} $2.4$' 0.024, '\fontsize{8}{8} $2.8$' 0.028, '\fontsize{8}{8} $3.2$' 0.032, '\fontsize{8}{8} $3.6$' 0.036) offset 0.6,-0.1

#set xtics ('\fontsize{8}{8} $0.06$' 0.06,'\fontsize{8}{8} $0.07$' 0.07, '\fontsize{8}{8} $0.08$' 0.08, '\fontsize{8}{8} $1.8$' 1.8, '\fontsize{8}{8} $2.0$' 2.0,'\fontsize{8}{8} $6.5$' 6.5) offset 0,0.25

set ytics ('\fontsize{8}{8} $1.0$' 0.010, '\fontsize{8}{8} $2.0$' 0.020, '\fontsize{8}{8} $3.0$' 0.030, '\fontsize{8}{8} $4.0$' 0.040, '\fontsize{8}{8} $5.0$' 0.050, '\fontsize{8}{8} $6.0$' 0.060, '\fontsize{8}{8} $7.0$' 0.070, '\fontsize{8}{8} $8.0$' 0.080) offset 0.75,-0.1

set xtics ('\fontsize{8}{8} $0.00$' 0.00,'\fontsize{8}{8} $0.02$' 0.02, '\fontsize{8}{8} $0.04$' 0.04, '\fontsize{8}{8} $0.06$' 0.06, '\fontsize{8}{8} $0.08$' 0.08,'\fontsize{8}{8} $0.10$' 0.10) offset -0.5,0.45

plot f(x) w l lt 1 lw 3 dt 3 lc rgb 'black' ti 'fit', "../../data/SpinSplitting_vs_Lx.txt" u (1.0)/($1):(2.0*$2) w lp lw 2.5 lt 1 pt 1 ps 1.5 lc rgb '#2e86c1' ti 'data'


#############################################################
unset multiplot
