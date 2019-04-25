gcc -g -Wall -Wformat-truncation=0 -O3 main.c mesh.c effectiveparticlepotential.c meshbiasfuncs.c utilities.c initialize.c parseparameters.c errorcheck.c -lm -lgsl -lgslcblas -o ../LatticeMembrane;
#  -no-pie -pg
#-g

