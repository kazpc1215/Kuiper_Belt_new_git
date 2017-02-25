reset
set term aqua dashed

set xl "semi-major axis[AU]" font "Times-Roman,30"
set yl "N" font "Times-Roman,30"
set xtics font "Times-Roman,30"
set ytics font "Times-Roman,30"
set title font "Times-Roman,25"
#set key font "Times-Roman,25"
set xr [20:70]
set yr [0:40]
unset key
set lmargin 11
set tmargin 2
set title offset 0,-1.5
set yl offset -1,0

n = 8

NEPfile = "./N100_1/Neptune.dat"
NEP_AXIS = system("cat " . NEPfile . " | awk \'NR==" . n . "{printf(\"%f\", $3" . ")}\'")
#NEP_AXIS = system("cat " . NEPfile . " | awk \'END{printf(\"%f\", $3" . ")}\'")
TIME = system("cat " . NEPfile . " | awk \'NR==" . n . "{printf(\"%e\", $1" . ")}\'")
#TIME = system("cat " . NEPfile . " | awk \'END{printf(\"%e\", $1" . ")}\'")
TIME = TIME/2.0/pi/1000000.0 


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



NEP_RESO21 = 2.0**(2.0/3.0)*NEP_AXIS
NEP_RESO32 = (3.0/2.0)**(2.0/3.0)*NEP_AXIS
NEP_RESO43 = (4.0/3.0)**(2.0/3.0)*NEP_AXIS
NEP_RESO54 = (5.0/4.0)**(2.0/3.0)*NEP_AXIS




set arrow 1 from NEP_RESO21,0 to NEP_RESO21,25 nohead dt 2 lw 2 lc rgb "gray50"
set arrow 2 from NEP_RESO32,0 to NEP_RESO32,25 nohead dt 2 lw 2 lc rgb "gray50"
set arrow 3 from NEP_RESO43,0 to NEP_RESO43,25 nohead dt 2 lw 2 lc rgb "gray50"
set arrow 4 from NEP_RESO54,0 to NEP_RESO54,25 nohead dt 2 lw 2 lc rgb "gray50"

set arrow 5 from NU5,0 to NU5,30 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 6 from NU6,0 to NU6,30 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 7 from NU7,0 to NU7,30 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 8 from NU8,0 to NU8,30 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 9 from NU15,0 to NU15,35 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 10 from NU16,0 to NU16,35 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 11 from NU17,0 to NU17,35 nohead dt 4 lw 2 lc rgb "gray50"
set arrow 12 from NU18,0 to NU18,35 nohead dt 4 lw 2 lc rgb "gray50"




set label "2:1" center at NEP_RESO21,25.5 font "Times-Roman,25"
set label "3:2" center at NEP_RESO32,27 font "Times-Roman,25"
set label "4:3" center at NEP_RESO43,25.5 font "Times-Roman,25"
set label "5:4" center at NEP_RESO54,27 font "Times-Roman,25"

set label "{/Symbol n}_5" center at NU5,31 font "Times-Roman,25"
set label "{/Symbol n}_6" center at NU6,32.5 font "Times-Roman,25"
set label "{/Symbol n}_7" center at NU7,32.5 font "Times-Roman,25"
set label "{/Symbol n}_8" center at NU8,31 font "Times-Roman,25"
set label "{/Symbol n}_{15}" center at NU15,36 font "Times-Roman,25"
set label "{/Symbol n}_{16}" center at NU16,36 font "Times-Roman,25"
set label "{/Symbol n}_{17}" center at NU17,36 font "Times-Roman,25"
set label "{/Symbol n}_{18}" center at NU18,36 font "Times-Roman,25"

if (n == 8){
set title sprintf("Time %.3f [Myr], Total 368",TIME)
plot "kuiperbelt_histogram_20Myr.dat" w boxes
}
if (n == 7){
set title sprintf("Time %.3f [Myr], Total 454",TIME)
plot "kuiperbelt_histogram_10Myr.dat" w boxes
}
if (n == 6){
set title sprintf("Time %.3f [Myr], Total 882",TIME)
plot "kuiperbelt_histogram_1Myr.dat" w boxes
}
if (n == 5){
set title sprintf("Time %.3f [Myr], Total 1000",TIME)
plot "kuiperbelt_histogram_100kyr.dat" w boxes
}
if (n == 4){
set title sprintf("Time %.3f [Myr], Total 1000",TIME)
plot "kuiperbelt_histogram_10kyr.dat" w boxes
}
if (n == 3){
set title sprintf("Time %.3f [Myr], Total 1000",TIME)
plot "kuiperbelt_histogram_1kyr.dat" w boxes
}
if (n == 2){
set title sprintf("Time %.3f [Myr], Total 1000",TIME)
plot "kuiperbelt_histogram_0yr.dat" w boxes
}

