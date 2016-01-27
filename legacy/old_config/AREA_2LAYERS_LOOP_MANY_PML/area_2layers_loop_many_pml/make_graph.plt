#! /usr/bin/gnuplot -persist
#set terminal png size 480, 360 enhanced font ',8'
#set terminal png size 600, 400 enhanced font ',8'
set terminal png size 800, 500 enhanced
set grid
set key top right
#set logscale y
set xlabel "x"

set style line 1 lt 1 pt 0 lw 2
set style line 2 lt 2 pt 0 lw 2
set style line 3 lt 3 pt 0 lw 2
set style line 4 lt 4 pt 0 lw 2
set style line 5 lt 5 pt 0 lw 2
set style line 6 lt 6 pt 0 lw 2
set style line 7 lt 7 pt 0 lw 2
set style line 8 lt 8 pt 0 lw 2
set style line 9 lt 9 pt 0 lw 2
set style line 10 lt 10 pt 0 lw 2

# ExR_10 2
# ExI_10 3
# EyR_10 4
# EyI_10 5
# EzR_10 6
# EzI_10 7
# ExR_-10 8
# ExI_-10 9
# EyR_-10 10
# EyI_-10 11
# EzR_-10 12
# EzI_-10 13
# smooth sbezier

set ylabel "EyR, z=10, y=0"
set output "EyR_10.png"
plot "line_pml_500.txt" using 1:4 title "PML в 500м от центра" with linespoints linestyle 2, \
     "line_pml_600.txt" using 1:4 title "PML в 600м от центра" with linespoints linestyle 3, \
     "line_pml_700.txt" using 1:4 title "PML в 700м от центра" with linespoints linestyle 4, \
     "line_pml_800.txt" using 1:4 title "PML в 800м от центра" with linespoints linestyle 5, \
     "line_pml_900.txt" using 1:4 title "PML в 900м от центра" with linespoints linestyle 6, \
     "line_pml_1000.txt" using 1:4 title "PML в 1000м от центра" with linespoints linestyle 7, \
     "line_pml_1100.txt" using 1:4 title "PML в 1100м от центра" with linespoints linestyle 8, \
     "line_pml_1200.txt" using 1:4 title "PML в 1200м от центра" with linespoints linestyle 9, \
     "line_pml_1300.txt" using 1:4 title "PML в 1300м от центра" with linespoints linestyle 10, \
     "line_std_600.txt" using 1:4 title "Большой бак" with linespoints linestyle 1, \


set ylabel "EyI, z=10, y=0"
set output "EyI_10.png"
plot "line_pml_500.txt" using 1:5 title "PML в 500м от центра" with linespoints linestyle 2, \
     "line_pml_600.txt" using 1:5 title "PML в 600м от центра" with linespoints linestyle 3, \
     "line_pml_700.txt" using 1:5 title "PML в 700м от центра" with linespoints linestyle 4, \
     "line_pml_800.txt" using 1:5 title "PML в 800м от центра" with linespoints linestyle 5, \
     "line_pml_900.txt" using 1:5 title "PML в 900м от центра" with linespoints linestyle 6, \
     "line_pml_1000.txt" using 1:5 title "PML в 1000м от центра" with linespoints linestyle 7, \
     "line_pml_1100.txt" using 1:5 title "PML в 1100м от центра" with linespoints linestyle 8, \
     "line_pml_1200.txt" using 1:5 title "PML в 1200м от центра" with linespoints linestyle 9, \
     "line_pml_1300.txt" using 1:5 title "PML в 1300м от центра" with linespoints linestyle 10, \
     "line_std_600.txt" using 1:5 title "Большой бак" with linespoints linestyle 1, \


set ylabel "EyR, z=-10, y=0"
set output "EyR_-10.png"
plot "line_pml_500.txt" using 1:10 title "PML в 500м от центра" with linespoints linestyle 2, \
     "line_pml_600.txt" using 1:10 title "PML в 600м от центра" with linespoints linestyle 3, \
     "line_pml_700.txt" using 1:10 title "PML в 700м от центра" with linespoints linestyle 4, \
     "line_pml_800.txt" using 1:10 title "PML в 800м от центра" with linespoints linestyle 5, \
     "line_pml_900.txt" using 1:10 title "PML в 900м от центра" with linespoints linestyle 6, \
     "line_pml_1000.txt" using 1:10 title "PML в 1000м от центра" with linespoints linestyle 7, \
     "line_pml_1100.txt" using 1:10 title "PML в 1100м от центра" with linespoints linestyle 8, \
     "line_pml_1200.txt" using 1:10 title "PML в 1200м от центра" with linespoints linestyle 9, \
     "line_pml_1300.txt" using 1:10 title "PML в 1300м от центра" with linespoints linestyle 10, \
     "line_std_600.txt" using 1:10 title "Большой бак" with linespoints linestyle 1, \


set ylabel "EyI, z=-10, y=0"
set output "EyI_-10.png"
plot "line_pml_500.txt" using 1:11 title "PML в 500м от центра" with linespoints linestyle 2, \
     "line_pml_600.txt" using 1:11 title "PML в 600м от центра" with linespoints linestyle 3, \
     "line_pml_700.txt" using 1:11 title "PML в 700м от центра" with linespoints linestyle 4, \
     "line_pml_800.txt" using 1:11 title "PML в 800м от центра" with linespoints linestyle 5, \
     "line_pml_900.txt" using 1:11 title "PML в 900м от центра" with linespoints linestyle 6, \
     "line_pml_1000.txt" using 1:11 title "PML в 1000м от центра" with linespoints linestyle 7, \
     "line_pml_1100.txt" using 1:11 title "PML в 1100м от центра" with linespoints linestyle 8, \
     "line_pml_1200.txt" using 1:11 title "PML в 1200м от центра" with linespoints linestyle 9, \
     "line_pml_1300.txt" using 1:11 title "PML в 1300м от центра" with linespoints linestyle 10, \
     "line_std_600.txt" using 1:11 title "Большой бак" with linespoints linestyle 1, \

