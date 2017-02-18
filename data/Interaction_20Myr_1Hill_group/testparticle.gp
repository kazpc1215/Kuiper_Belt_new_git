reset
set term aqua dashed

set xl "semi-major axis[AU]" font "Times-Roman,30"
set yl "ecc" font "Times-Roman,30"
set xtics font "Times-Roman,30"
set ytics font "Times-Roman,30"
set title font "Times-Roman,25"
#set key font "Times-Roman,25"
unset key
set lmargin 11
set tmargin 2
set title offset 0,-1.5
set yl offset -1,0

n = 2
j = 0
while (n<=2){
unset label
unset arrow

NUfile = "nu_NEP.dat"
NU8 = system("cat " . NUfile . " | awk \'NR==" . n . "{printf(\"%f\", $2" . ")}\'")
NU18 = system("cat " . NUfile . " | awk \'NR==" . n . "{printf(\"%f\", $3" . ")}\'")
NUfile = "nu_JUP.dat"
NU5 = system("cat " . NUfile . " | awk \'NR==" . n . "{printf(\"%f\", $2" . ")}\'")
NU15 = system("cat " . NUfile . " | awk \'NR==" . n . "{printf(\"%f\", $3" . ")}\'")
NUfile = "nu_SAT.dat"
NU6 = system("cat " . NUfile . " | awk \'NR==" . n . "{printf(\"%f\", $2" . ")}\'")
NU16 = system("cat " . NUfile . " | awk \'NR==" . n . "{printf(\"%f\", $3" . ")}\'")
NUfile = "nu_URA.dat"
NU7 = system("cat " . NUfile . " | awk \'NR==" . n . "{printf(\"%f\", $2" . ")}\'")
NU17 = system("cat " . NUfile . " | awk \'NR==" . n . "{printf(\"%f\", $3" . ")}\'")


NEPfile = "./N100_1/Neptune.dat"
NEP_AXIS = system("cat " . NEPfile . " | awk \'NR==" . n . "{printf(\"%f\", $3" . ")}\'")
#NEP_AXIS = system("cat " . NEPfile . " | awk \'END{printf(\"%f\", $3" . ")}\'")
NEP_RH = system("cat " . NEPfile . " | awk \'NR==" . n . "{printf(\"%f\", $8" . ")}\'")
#NEP_RH = system("cat " . NEPfile . " | awk \'END{printf(\"%f\", $8" . ")}\'")
TIME = system("cat " . NEPfile . " | awk \'NR==" . n . "{printf(\"%e\", $1" . ")}\'")
#TIME = system("cat " . NEPfile . " | awk \'END{printf(\"%e\", $1" . ")}\'")
TIME = TIME/2.0/pi/1000000.0 
set title sprintf("Time %.3f [Myr]",TIME)

unset surface
set contour
set cntrparam bspline
set cntrparam order 10
#set cntrparam levels incremental 2.8, 0.05, 3.2
set samples 1000
set isosamples 100
unset surface
set view 0,0
set xr [0:70]
set yr [0:0.5]
set zr [0:5]
I = 0.0
#
NEP_cont1RH = NEP_AXIS/(NEP_AXIS+NEP_RH) + 2.0*sqrt((NEP_AXIS+NEP_RH)/NEP_AXIS)
NEP_cont4RH = NEP_AXIS/(NEP_AXIS+4.0*NEP_RH) + 2.0*sqrt((NEP_AXIS+4.0*NEP_RH)/NEP_AXIS)
NEP_cont10RH = NEP_AXIS/(NEP_AXIS+10.0*NEP_RH) + 2.0*sqrt((NEP_AXIS+10.0*NEP_RH)/NEP_AXIS)
set cntrparam levels discrete NEP_cont1RH, NEP_cont4RH, NEP_cont10RH
set table "Jacobi_contour_NEP" 
splot NEP_AXIS/x+2.0*sqrt(x/NEP_AXIS)*sqrt(1.0-y*y)*cos(I)
unset table

NEP_RESO21 = 2.0**(2.0/3.0)*NEP_AXIS
NEP_RESO32 = (3.0/2.0)**(2.0/3.0)*NEP_AXIS
NEP_RESO43 = (4.0/3.0)**(2.0/3.0)*NEP_AXIS
NEP_RESO54 = (5.0/4.0)**(2.0/3.0)*NEP_AXIS

set arrow 1 from NEP_RESO21,0 to NEP_RESO21,0.4 nohead dt 2 lw 2 lc rgb "gray50"
set arrow 2 from NEP_RESO32,0 to NEP_RESO32,0.4 nohead dt 2 lw 2 lc rgb "gray50"
set arrow 3 from NEP_RESO43,0 to NEP_RESO43,0.4 nohead dt 2 lw 2 lc rgb "gray50"
set arrow 4 from NEP_RESO54,0 to NEP_RESO54,0.4 nohead dt 2 lw 2 lc rgb "gray50"

#set arrow 5 from NU5,0 to NU5,0.45 nohead dt 4 lw 2 lc rgb "gray50"
#set arrow 6 from NU6,0 to NU6,0.45 nohead dt 4 lw 2 lc rgb "gray50"
#set arrow 7 from NU7,0 to NU7,0.45 nohead dt 4 lw 2 lc rgb "gray50"
#set arrow 8 from NU8,0 to NU8,0.45 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 9 from NU15,0 to NU15,0.45 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 10 from NU16,0 to NU16,0.45 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 11 from NU17,0 to NU17,0.45 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 12 from NU18,0 to NU18,0.45 nohead dt 4 lw 2 lc rgb "gray50"


set label "2:1" center at NEP_RESO21,0.41 font "Times-Roman,25"
set label "3:2" center at NEP_RESO32,0.43 font "Times-Roman,25"
set label "4:3" center at NEP_RESO43,0.41 font "Times-Roman,25"
set label "5:4" center at NEP_RESO54,0.43 font "Times-Roman,25"

#set label "{/Symbol n}_5" center at NU5,0.47 font "Times-Roman,25"
#set label "{/Symbol n}_6" center at NU6,0.49 font "Times-Roman,25"
#set label "{/Symbol n}_7" center at NU7,0.49 font "Times-Roman,25"
#set label "{/Symbol n}_8" center at NU8,0.47 font "Times-Roman,25"
set label "{/Symbol n}_{15}" center at NU15,0.47 font "Times-Roman,25"
set label "{/Symbol n}_{16}" center at NU16,0.47 font "Times-Roman,25"
set label "{/Symbol n}_{17}" center at NU17,0.47 font "Times-Roman,25"
set label "{/Symbol n}_{18}" center at NU18,0.47 font "Times-Roman,25"

set xr [20:70]
set yr [0:0.5]

#plot for [i=1:100] sprintf("./N100_1/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_2/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_3/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_4/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_5/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_6/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_7/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_8/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_9/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_10/test_particle%03d.dat",i) every ::j::j u 3:2 w p pt 5 lc rgb "black",\
"./N100_1/Jupiter.dat" every ::j::j u 3:2:8 w circle lt 1,\
"./N100_1/Saturn.dat" every ::j::j u 3:2:8 w circle lt 2,\
"./N100_1/Uranus.dat" every ::j::j u 3:2:8 w circle lt 3,\
"./N100_1/Neptune.dat" every ::j::j u 3:2:8 w circle lt 4,\
1.0 - NEP_AXIS/x lt 4 dt 2
#"Jacobi_contour_NEP" u 1:2 w l dt 2 lt 4 lw 2

set yl "inc[rad]"
plot for [i=1:100] sprintf("./N100_1/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_2/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_3/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_4/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_5/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_6/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_7/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_8/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_9/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
for [i=1:100] sprintf("./N100_10/test_particle%03d.dat",i) every ::j::j u 3:5 w p pt 5 lc rgb "black",\
"./N100_1/Jupiter.dat" every ::j::j u 3:5:8 w circle lt 1,\
"./N100_1/Saturn.dat" every ::j::j u 3:5:8 w circle lt 2,\
"./N100_1/Uranus.dat" every ::j::j u 3:5:8 w circle lt 3,\
"./N100_1/Neptune.dat" every ::j::j u 3:5:8 w circle lt 4

n = n + 1
j = j + 1
}



