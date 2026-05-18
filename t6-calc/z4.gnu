set terminal pngcairo size 700,700 enhanced font 'Arial,11'
set output 'z4.png'
set title 'Задание 4: полярная фигура,N=10000' font 'Arial,12'
set xlabel 'x'
set ylabel 'y'
set grid lw 1 lc rgb '#cccccc'
set size ratio 1
set key top right
plot 'z4_out.dat' using 1:2 with points pt 1 ps 0.3 lc rgb '#AAAAAA' title 'вне',\
     'z4_in.dat'  using 1:2 with points pt 1 ps 0.3 lc rgb '#CC2200' title 'внутри',\
     'z4_crv.dat' using 1:2 with lines lw 2 lc rgb '#2266CC' title 'ρ(φ)'
