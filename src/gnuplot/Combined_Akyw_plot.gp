set terminal epslatex size 4,3.5 standalone color colortext 8
set output 'OS_resolved_Akyw_for_Vo_9p5_36x36_SSW.tex'

set border linewidth 2.75 lc rgb 'white'
set tics textcolor rgb 'black'

#########################################################
set xrange [-pi:pi]
set yrange [-0.6:0.6]

set xzeroaxis lw 4 lt 1 dt 3 lc rgb 'white'
set ytics nomirror
set xtics nomirror

set pm3d map
set pm3d interpolate 0,0
unset colorbox

set multiplot
#########################################################

set lmargin at screen 0.08
set rmargin at screen 0.52
set bmargin at screen 0.545
set tmargin at screen 0.98

set ylabel rotate by 0 '\scalebox{1.4}{$\frac{\omega}{m}$}' offset 6.1,0

set ytics ('\fontsize{9}{9} $-0.2$' -0.2, '\fontsize{9}{9} $-0.4$' -0.4, '\fontsize{9}{9} $-0.6$' -0.6, '\fontsize{9}{9} $0.0$' 0.0, '\fontsize{9}{9} $0.2$' 0.2, '\fontsize{9}{9} $0.4$' 0.4, '\fontsize{9}{9} $0.6$' 0.6) offset 1.1,0 

set xtics ('' -pi, '' -0.5*pi, '' 0.0, '' 0.5*pi, '' pi) offset 0,0
set label 1 '\scalebox{1.25}{(a)}' at -3.1,0.53 tc rgb 'white' front
set label 2 '\scalebox{1.5}{$s,\uparrow$}' at 2.2,0.54 tc rgb 'white' front

splot "../../data/Akyw_s_for_Vo9.5_A0.3_B0.5_D0.0_HE0.3.txt" u 1:2:3 w pm3d ti ''

unset label 1
unset label 2
unset ylabel 
#########################################################

set lmargin at screen 0.55
set rmargin at screen 0.99
set bmargin at screen 0.545
set tmargin at screen 0.98

set ytics ('' -0.2, '' -0.4, '' -0.6, '' 0.0, '' 0.2, '' 0.4, '' 0.6) offset 0.75,0
set label 1 '\scalebox{1.25}{(b)}' at -3.1,0.53 tc rgb 'white' front
set label 2 '\scalebox{1.5}{$s,\downarrow$}' at 2.2,0.54 tc rgb 'white' front

splot "../../data/Akyw_s_for_Vo9.5_A0.3_B0.5_D0.0_HE0.3.txt" u 1:2:4 w pm3d ti ''

unset label 1
unset label 2
#########################################################

set lmargin at screen 0.08
set rmargin at screen 0.52
set bmargin at screen 0.08
set tmargin at screen 0.515

set ylabel rotate by 0 '\scalebox{1.4}{$\frac{\omega}{m}$}' offset 6.1,-0.1
set xlabel '\fontsize{12}{12} $k_y a$' offset 0.5,1.1

set ytics ('\fontsize{9}{9} $-0.2$' -0.2, '\fontsize{9}{9} $-0.4$' -0.4, '\fontsize{9}{9} $-0.6$' -0.6, '\fontsize{9}{9} $0.0$' 0.0, '\fontsize{9}{9} $0.2$' 0.2, '\fontsize{9}{9} $0.4$' 0.4, '\fontsize{9}{9} $0.6$' 0.6) offset 1.1,-0.2
set xtics ('\fontsize{9}{9} $-\pi$' -1.0*pi, '\fontsize{9}{9} $-\pi/2$' -0.5*pi, '\fontsize{9}{9} $0$' 0.0, '\fontsize{9}{9} $\pi/2$' 0.5*pi, '\fontsize{9}{9} $\pi$' 1.0*pi) offset -0.3,0.7

set label 1 '\scalebox{1.25}{(c)}' at -3.1,0.53 tc rgb 'white' front
set label 2 '\scalebox{1.5}{$p,\uparrow$}' at 2.2,0.54 tc rgb 'white' front

splot "../../data/Akyw_p_for_Vo9.5_A0.3_B0.5_D0.0_HE0.3.txt" u 1:2:3 w pm3d ti ''

unset label 1
unset label 2
unset ylabel
#########################################################

set lmargin at screen 0.55
set rmargin at screen 0.99
set bmargin at screen 0.08
set tmargin at screen 0.515

set ytics ('' -0.2, '' -0.4, '' -0.6, '' 0.0, '' 0.2, '' 0.4, '' 0.6) offset 0.75,0
set label 1 '\scalebox{1.25}{(d)}' at -3.1,0.53 tc rgb 'white' front
set label 2 '\scalebox{1.5}{$p,\downarrow$}' at 2.2,0.54 tc rgb 'white' front

splot "../../data/Akyw_p_for_Vo9.5_A0.3_B0.5_D0.0_HE0.3.txt" u 1:2:4 w pm3d ti ''

unset label 1
unset label 2



unset multiplot
