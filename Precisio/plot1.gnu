#Script histograma
reset
set term png
#Nom del fitxer de sortida
set output "res_final.png"

#Titols de la gràfica i eixos
set title "E pot"
set xlabel "T"
set ylabel "E pot"

#Movem informació dades
#set key 

#Escala logaritmica en un eix(o dos)
#set logscale x

#Cal·librem els eixos
set xrange[:]
set yrange[:]
#t=1.0
#e=1.0
#f(x)=((1*exp(-x*1))/(sqrt(1-exp(-x))*(1-exp(-1*x))))-(1/(x**(3/2)))
plot 'result.dat' i 0 u 1:2 w p,'' i 1 u 1:2 w p

