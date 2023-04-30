reset
set term gif size 1000,600 animate delay 10 loop 0 optimize
set output "funcio_ona.gif"

# leeremos el numero de bloques de manera automatica

datafile ="program-res1.dat"
stats datafile
numblo=STATS_blocks

# pone el titulo
set title "T(x,t)"

# fijamos los rangos del canvas, para que no se reajusten
set xrange[-17:17]
set yrange[0:1]

#fijamos los títulos
set xlabel "x(m)"
set ylabel "T(x,t) (ºC)"

#bucle
do for[i=0:numblo-2]{

#escribe la etiqueta del tiempo
set label 2 sprintf('Time: %9.3f (sec)',i*0.01) at 10,160 right front font 'Verdana,12'

#dibuja
plot datafile index i u 2:3 w l t"T(x,t)","" index 0 u 2:3 w l t"T_0(x)", "potential.dat" u 1:2 w l t'V(x)'

#borra la etiqueta
unset label 2
}
reset
set term png
set output"transmission_coeficient.png"
plot 'program-coeficients.dat' u 1:2 w l t'transmission', "" u 1:3 w l t"reflexion"
