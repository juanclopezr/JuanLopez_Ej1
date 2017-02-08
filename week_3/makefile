v.png: plotter.py output.txt
	python plotter.py
	rm -f output.txt placas

output.txt : a.out 
	mpiexec -n 4 ./a.out > output.txt

a.out : placa.c
	mpicc placa.c
