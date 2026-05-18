set terminal pngcairo size 900,600 enhanced font 'Arial,11'
set output 'z2.png'
set title 'Задание 2: Монте-Карло,интеграл,N=10000' font 'Arial,12'
set xlabel 'x'
set ylabel 'y'
set grid lw 1 lc rgb '#cccccc'
set key top right
plot 'z2_out.dat' using 1:2 with points pt 1 ps 0.3 lc rgb '#AAAAAA' title 'вне',\
     'z2_in.dat'  using 1:2 with points pt 1 ps 0.3 lc rgb '#CC2200' title 'внутри',\
     'z2_crv.dat' using 1:2 with lines lw 2 lc rgb '#2266CC' title 'sqrt(29-14cos²x)'
