set terminal pngcairo size 900,600 enhanced font 'Arial,11'
set output 'z1.png'
set title 'Задание 1: Монте-Карло,треугольник n=14\nS_MC=99.04 (точно=100),N=10000' font 'Arial,12'
set xlabel 'x'
set ylabel 'y'
set grid lw 1 lc rgb '#cccccc'
set xrange [-0.5:21]
set yrange [-0.5:11.5]
set key top right
plot 'z1_out.dat' using 1:2 with points pt 1 ps 0.3 lc rgb '#BBBBBB' title 'вне треугольника',\
     'z1_in.dat'  using 1:2 with points pt 1 ps 0.3 lc rgb '#CC2200' title 'внутри',\
     'z1_crv.dat' using 1:2 with lines lw 2.5 lc rgb '#2266CC' title 'f1(x)=10x/14',\
     'z1_crv.dat' using 1:3 with lines lw 2.5 lc rgb '#228833' title 'f2(x)=10(x-20)/(14-20)',\
     'z1_vx.dat' using 1:2 with points pt 7 ps 2 lc rgb '#000000' title 'вершины'
