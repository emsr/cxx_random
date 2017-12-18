gnuplot

set termoption enhanced

set xzeroaxis
set yzeroaxis
set grid

set title "Poisson sampling"
set xlabel "k"
plot [90.0:230.0][0.0:3250000.0] \
                    "test_poisson.txt" index 0 using 1:2 with lines title "poop"

set title "Poisson sampling"
set xlabel "k"
plot [150.0:165.0][2900000.0:3250000.0] \
                    "test_poisson.txt" index 0 using 1:2 with lines title "poop"

set title "Poisson sampling"
set xlabel "k"
plot [90.0:230.0][0.0:0.033] \
                    "test_poisson.txt" index 0 using 1:3 with lines title "poop"

set title "Poisson sampling"
set xlabel "k"
plot [150.0:165.0][0.025:0.033] \
                    "test_poisson.txt" index 0 using 1:3 with lines title "poop"
