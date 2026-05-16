set terminal pngcairo size 800,600
set output 'comp_num_analit'
set termoption enhanced

set title 'Solução numérica comparada à analítica da literatura'
set xlabel 'Tempo'
set xrange [0:11]
set ylabel 'Velocidade'
set yrange [0:1.1]
set key right bottom
set grid

plot\
    'sol_num01.dat' using 1:2 with lines lw 2 title 'Solução numérica com Res = 0.1', \
    'sol_num01.dat' using 1:3 with lines lw 2 title 'Solução analítica com Res = 0.1', \
    'sol_num04.dat' using 1:2 with lines lw 2 title 'Solução numérica com Res = 0.4', \
    'sol_num04.dat' using 1:3 with lines lw 2 title 'Solução analítica com Res = 0.4', \
    'sol_num07.dat' using 1:2 with lines lw 2 title 'Solução numérica com Res = 0.7', \
    'sol_num07.dat' using 1:3 with lines lw 2 title 'Solução analítica com Res = 0.7', \
    'sol_num1.dat' using 1:2 with lines lw 2 title 'Solução numérica com Res = 1.0', \
    'sol_num1.dat' using 1:3 with lines lw 2 title 'Solução analítica com Res = 1.0'
