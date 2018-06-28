default: heisenberg3d_4layer.c plot_h3d_4layer.py
	gcc -fPIC -shared -o heisenberg3d_4layer.so -g -I /usr/include -L /usr/include heisenberg3d_4layer.c -lm -lgsl -lgslcblas 
	nohup python3 plot_h3d_4layer.py &
