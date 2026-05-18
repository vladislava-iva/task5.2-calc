set terminal pngcairo size 700,700 enhanced font 'Arial,11'
set output 'z3.png'
set title 'Задание 3: π ≈ 3.1580,N=10000' font 'Arial,12'
set xlabel 'x'
set ylabel 'y'
set grid lw 1 lc rgb '#cccccc'
set size ratio 1
set xrange [-15.0000:15.0000]
set yrange [-15.0000:15.0000]
set key top right
plot 'z3_out.dat' using 1:2 with points pt 1 ps 0.3 lc rgb '#AAAAAA' title 'вне круга',\
     'z3_in.dat'  using 1:2 with points pt 1 ps 0.3 lc rgb '#CC2200' title 'в круге',\
     'z3_circ.dat' using 1:2 with lines lw 2 lc rgb '#2266CC' title 'R=14.0000'
