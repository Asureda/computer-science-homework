reset
set term png
#Nom del fitxer de sortida
set output "NN_cos.png"

#Titols de la gràfica i eixos
set title "Interpolation cos(x)"
set xlabel "x"
set ylabel "y"

#Movem informació dades
#set key 

#Escala logaritmica en un eix(o dos)
#set logscale y

#Cal·librem els eixos
set xrange[0.5:10]
set yrange[:]
plot "known.txt" i 0 u 1:2 w p t"cos(x)","" "interp.txt" u 1:2 w l t"interpolation", cos(x) w l t"cos(x)"