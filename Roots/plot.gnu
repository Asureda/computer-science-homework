#Script histograma
reset
set term png
#Nom del fitxer de sortida
set output "Absolute_Error.png"

#Titols de la gràfica i eixos
set title "Absolute Error"
set xlabel "n"
set ylabel "Abs Error"

#Movem informació dades
#set key 

#Escala logaritmica en un eix(o dos)
set logscale y

#Cal·librem els eixos
set xrange[:]
set yrange[:]
#a=0.5+sin(1)-sin(0)
plot "method_convergence.dat" i 0 u 1:2 w lp t"Bisection", "" i 1 u 1:2 w lp t"Secant","" i 2 u 1:2 w lp t"Regula Falsi",  "" i 3 u 1:2 w lp t"Newton-Raphson","" i 4 u 1:2 w lp t"Mullers"


#Script histograma
reset
set term png
#Nom del fitxer de sortida
set output "Relative_Error.png"

#Titols de la gràfica i eixos
set title "Relative Error"
set xlabel "n"
set ylabel "Rel Error"

#Movem informació dades
#set key 

#Escala logaritmica en un eix(o dos)
set logscale y

#Cal·librem els eixos
set xrange[:]
set yrange[:]
#a=0.5+sin(1)-sin(0)
plot "method_convergence.dat" i 0 u 1:3 w lp t"Bisection", "" i 1 u 1:3 w lp t"Secant","" i 2 u 1:3 w lp t"Regula Falsi",  "" i 3 u 1:3 w lp t"Newton-Raphson","" i 4 u 1:3 w lp t"Mullers"