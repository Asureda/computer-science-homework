#Script histograma
reset
set term png
#Nom del fitxer de sortida
set output "Integral.png"

#Titols de la gràfica i eixos
set title "E pot"
set xlabel "T"
set ylabel "E pot"

#Movem informació dades
#set key 

#Escala logaritmica en un eix(o dos)
#set logscale x

#Cal·librem els eixos
set xrange[0:1000]
set yrange[:]
#t=1.0
#e=1.0
#f(x)=((1*exp(-x*1))/(sqrt(1-exp(-x))*(1-exp(-1*x))))-(1/(x**(3/2)))
plot 'integral_conv.dat' i 1 u 2:1 w lp

#Script histograma
reset
set term png
#Nom del fitxer de sortida
set output "Funcio_integrand.png"

#Titols de la gràfica i eixos
set title "E pot"
set xlabel "T"
set ylabel "E pot"

#Movem informació dades
#set key 

#Escala logaritmica en un eix(o dos)
#set logscale y

#Cal·librem els eixos
set xrange[0:20]
#set yrange[1:100000]
E=0.1
etha=0.5
f(x)=((etha*exp(-E*x))/(sqrt(1-exp(-x))*(1-exp(-etha*x))))-(1/(x**(3/2)))
plot 'integral_conv.dat' i 0 u 1:2 w lp , f(x) w l