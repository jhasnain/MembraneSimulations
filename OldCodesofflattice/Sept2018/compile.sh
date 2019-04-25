gcc -Wall -Wformat-truncation=0 -O3 main.c mesh.c meshbiasfuncs.c particles.c gcmc.c utilities.c initialize.c parseparameters.c bondinteractions.c errorcheck.c -lm -lgsl -lgslcblas -o ../MembraneFusion;
# -no-pie -pg 
#-g

