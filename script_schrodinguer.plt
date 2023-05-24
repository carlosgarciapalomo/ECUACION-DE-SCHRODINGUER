#Script que representa la probabilidad de la función de onda a lo largo del tiempo

reset session

stats 'Probabilidad.txt' nooutput
Columnas=STATS_columns
Filas=STATS_records

izqd=2*(Columnas/2)/5
drch=3*(Columnas/2)/5

Delay=Filas/60

set yrange [0:1]
set xlabel 'Posición'
set ylabel 'Probabilidad'

unset key

set style fill solid border -1
set object rectangle from izqd,0 to drch,0.7 fc rgb "blue" back

set terminal gif size 800,800 animate delay Delay optimize
set output 'Probabilidad.gif'

do for [j=0:Filas-1]{
	set table $Extract
		plot for [i=1:Columnas/2] 'Probabilidad.txt' u (column(2*i-1)):(column(2*i)) every ::j::j w table
	unset table
	plot $Extract u 1:2 w boxes 
}


set output 