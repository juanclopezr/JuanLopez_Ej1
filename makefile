placas.pdf : grafica.py potential.dat fieldx.dat fieldy.dat
	python grafica.py

potential.dat fieldx.dat fieldy.dat: a.out
	mpiexec -n 8 ./a.out

a.out : placas.c
	mpicc placas.c

clean :
	rm -f potential.dat fieldx.dat fieldy.dat a.out
