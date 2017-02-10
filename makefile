placas.pdf : grafica.py potential.dat fieldx.dat fieldy.dat
	python grafica.py
	rm -f potential.dat fieldx.dat fieldy.dat a.out

#submit : a.out
#	qsub submit.sh

a.out : placas.c
	mpicc placas.c
