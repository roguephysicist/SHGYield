#    	G N U P L O T
#    	Version 5.0 patchlevel 6    last modified 2017-03-18
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2017
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
set terminal pdfcairo transparent enhanced fontscale 0.7 size 7.00in, 5.00in
set output 'si111-shgyield.pdf'

set key top left height 0.4 width 1 samplen -1
set grid
set border lw 2

set xrange [ 2.5 : 5 ] noreverse nowriteback
set format y '%.2f'

LEFT = 0.13
TOP = 0.90
DIV = 0.40
SEP = 0.04

set multiplot layout 2,2 title 'SSHG Yield for the Si(111)(1x1):H surface'

unset xlabel
set format x ''
set ylabel '{R (10^{-20} x cm^{2}/W)}'
set yrange [ 0 : 3 ] noreverse nowriteback
set ytics 1
set tmargin at screen TOP
set bmargin at screen TOP-DIV+SEP
set lmargin at screen LEFT
set rmargin at screen LEFT+DIV-SEP
p 'si111-reference.out' u (2*$1):2 t 'R_{pP}' lw 2 lc rgb "#2aa198" w l

unset xlabel
set format x ''
unset ylabel
set yrange [ 0 : 0.4 ] noreverse nowriteback
set ytics 0.1
set tmargin at screen TOP
set bmargin at screen TOP-DIV+SEP
set lmargin at screen LEFT+DIV+SEP
set rmargin at screen LEFT+2*DIV
p 'si111-reference.out' u (2*$1):3 t 'R_{pS}' lw 2 lc rgb "#dc322f" w l

set xlabel "Two-photon energy (eV)" 
set format x '%.1f'
set ylabel '{R (10^{-20} x cm^{2}/W)}'
set yrange [ 0 : 0.1 ] noreverse nowriteback
set ytics 0.02
set tmargin at screen TOP-DIV-SEP
set bmargin at screen TOP-2*DIV
set lmargin at screen LEFT
set rmargin at screen LEFT+DIV-SEP
p 'si111-reference.out' u (2*$1):4 t 'R_{sP}' lw 2 lc rgb "#6c71c4" w l

set xlabel "Two-photon energy (eV)" 
set format x '%.1f'
unset ylabel
set yrange [ 0 : 0.1 ] noreverse nowriteback
set ytics 0.02
set tmargin at screen TOP-DIV-SEP
set bmargin at screen TOP-2*DIV
set lmargin at screen LEFT+DIV+SEP
set rmargin at screen LEFT+2*DIV
p 'si111-reference.out' u (2*$1):5 t 'R_{sS}' lw 2 lc rgb "#859900" w l

unset multiplot

#    EOF
