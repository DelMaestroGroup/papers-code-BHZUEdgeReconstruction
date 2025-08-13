set terminal epslatex size 3.5,1.75 standalone color colortext 8
set output 'Akyw_for_diff_eps_Vo_10p0_36x36_SSW.tex'

set border linewidth 2.75 lc rgb 'white'
set tics textcolor rgb 'black'

#########################################################
set xrange [-pi:pi]
set yrange [-0.6:0.6]

set xzeroaxis lw 3 lt 1 dt 3 lc rgb 'white'
set ytics nomirror
set xtics nomirror

set xlabel '\scalebox{1.1}{$k_y a$}' offset 0.5,1.4
set xtics ('\fontsize{7}{7} $-\pi$' -1.0*pi, '\fontsize{7}{7} $-\pi/2$' -0.5*pi, '\fontsize{7}{7} $0$' 0.0, '\fontsize{7}{7} $\pi/2$' 0.5*pi, '\fontsize{7}{7} $\pi$' 1.0*pi) offset -0.4,0.8

set pm3d map
set pm3d interpolate 0,0
unset colorbox

set multiplot
#########################################################

set lmargin at screen 0.08
set rmargin at screen 0.515
set bmargin at screen 0.13
set tmargin at screen 0.97

set ylabel rotate by 0 '\scalebox{1.1}{$\frac{\omega}{m}$}' offset 6.5,0

set ytics ('\fontsize{7}{7} $-0.2$' -0.2, '\fontsize{7}{7} $-0.4$' -0.4, '\fontsize{7}{7} $-0.6$' -0.6, '\fontsize{7}{7} $0.0$' 0.0, '\fontsize{7}{7} $0.2$' 0.2, '\fontsize{7}{7} $0.4$' 0.4, '\fontsize{7}{7} $0.6$' 0.6) offset 1.1,0 

set label 1 '\scalebox{1.1}{(a)}' at -3.1,0.53 tc rgb 'white' front
#set label 2 '\scalebox{1.1}{$\varepsilon_h=0.0$}' at 1.1,0.53 tc rgb 'white' front

splot "../data/Akyw_for_Vo10.0_A0.3_B0.5_D0.0_HE0.0.txt" u 1:2:($5*10) w pm3d ti ''

unset label 1
unset label 2
unset ylabel 
#########################################################

set lmargin at screen 0.555
set rmargin at screen 0.99
set bmargin at screen 0.13
set tmargin at screen 0.97

set xtics ('\fontsize{7}{7} $-\pi$' -1.0*pi, '\fontsize{7}{7} $-\pi/2$' -0.5*pi, '\fontsize{7}{7} $0$' 0.0, '\fontsize{7}{7} $\pi/2$' 0.5*pi, '\fontsize{7}{7} $\pi$' 1.0*pi) offset -0.2,0.8

set ytics ('' -0.2, '' -0.4, '' -0.6, '' 0.0, '' 0.2, '' 0.4, '' 0.6) offset 0.75,0
set label 1 '\scalebox{1.1}{(b)}' at -3.1,0.53 tc rgb 'white' front
#set label 2 '\scalebox{1.1}{$\varepsilon_h=0.3$}' at 1.1,0.53 tc rgb 'white' front

splot "../data/Akyw_for_Vo10.0_A0.3_B0.5_D0.0_HE0.3.txt" u 1:2:5 w pm3d ti ''

unset label 1
unset label 2
#########################################################

unset multiplot
