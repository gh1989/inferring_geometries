set term gif animate delay 3 size 400, 400
set output "trajectory.gif"
f(x,y) = cos(2*pi*x) + sin(2*pi*x) + cos(2*pi*y) + sin(2*pi*y)
g(x,y)=0
set palette rgbformulae 33,13,10
set sample 128
set pm3d map
set isosample 128
set palette defined (-1 "white", 0 "white", 1 "red")
do for [n=1:250] {
splot [-1:1][-1:1][-1:1] f(x,y) w pm3d,\
"data.txt" u 2:3:(g($2,$3)) every :::::n w lp t sprintf("n=%i", n) pt 0,\
"data.txt" u 4:5:(g($4,$5)) every :::::n w lp t sprintf("") pt 0,\
"data.txt" u 6:7:(g($6,$7)) every :::::n w lp t sprintf("") pt 0,\
"data.txt" u 8:9:(g($8,$9)) every :::::n w lp t sprintf("") pt 0,\
"data.txt" u 10:11:(g($10,$11)) every :::::n w lp t sprintf("") pt 0,\
"data.txt" u 12:13:(g($12,$13)) every :::::n w lp t sprintf("") pt 0,\
"data.txt" u 14:15:(g($14,$15)) every :::::n w lp t sprintf("") pt 0,\
"data.txt" u 16:17:(g($16,$17)) every :::::n w lp t sprintf("") pt 0,\
"data.txt" u 18:19:(g($18,$19)) every :::::n w lp t sprintf("") pt 0,\
"data.txt" u 20:21:(g($20,$21)) every :::::n w lp t sprintf("") pt 0\
}