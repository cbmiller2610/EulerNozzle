gcc -fPIC -Wall -shared -L /usr/local/lib -lgsl -lgslcblas -lm heat_solver.c -o heat_solver.so
