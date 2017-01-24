graph.pdf : graph.py random.dat
	python ghraph.py

random.dat : radioactive.x
	./radioactive.x

radioactive.x : radioactive.c
	gcc radioactive.c -o radioactive.x -lm
