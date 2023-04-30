#Script histograma
reset
set term png
#Nom del fitxer de sortida
set output "en_conv.png"

#Titols de la gràfica i eixos
set title "Convervence of e_1 with n max"
set xlabel "n"
set ylabel "E_1"

#Movem informació dades
#set key 

#Escala logaritmica en un eix(o dos)
#set logscale x

#Cal·librem els eixos
#set xrange[0:1000]
set yrange[0.3675:0.368]

plot 'en_conv.dat' u 1:2 w lp t'E_1'

#Script histograma
reset
set term png
#Nom del fitxer de sortida
set output "e_conv.png"

#Titols de la gràfica i eixos
set title "Convergence of e using E_1"
set xlabel "n"
set ylabel "e"

#Movem informació dades
#set key 

#Escala logaritmica en un eix(o dos)
#set logscale y

#Cal·librem els eixos
set yrange[2.71:2.72]
#set yrange[1:100000]

plot 'e_conv.dat' u 1:2 w lp t'e'