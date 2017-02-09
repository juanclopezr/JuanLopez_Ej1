#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

int L = 5, l = 2, d = 1, V0 = 100, m = 256, size, rank, m_y, min_range, max_range;

int transformer(int i, int j);
double *init(int x0, int x1, int y0, int y1, double *array);

int main(void)
{
    /*
     * Requirements 
     * */
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Request send_request, recv_request;
    MPI_Status status;
    
	int up, down, left, right, j=1, k = 0, n = 0;
    int x0, x1, y0, y1, N, i;
    double average, h = (double)L/m;
	
    m_y = m/size;
    N = 1000;//2*m*m_y;
    max_range = m_y*(rank + 1);
    min_range = m_y*(rank);
        
	x0 = m/2 - l/(h*2) - 1;
	x1 = m/2 + l/(h*2) - 1;        
	y0 = m/2 - d/(h*2) - 1;
	y1 = m/2 + d/(h*2) - 1;
    
    double *V = malloc((m*m_y+1)*sizeof(double));
    double *front_1 = malloc(m*sizeof(double));
    double *front_2 = malloc(m*sizeof(double));
    double *buff_1 = malloc(m*sizeof(double));
    double *buff_2 = malloc(m*sizeof(double));
    V = init(x0, x1, y0, y1, V);

    while (n < N)
    {		
        for(i=m_y*rank; i < (rank+1)*m_y; i++)
        {
            for(j=1;j < m; j++)
            {
                up = transformer(i-1, j);
                down = transformer(i+1, j);
                left = transformer(i, j-1);
                right = transformer(i, j+1);
		if(rank==0)
		  {
		    if(i==0)
		      {
			average = 0;
			V[transformer(i,j)] = average;
		      }
		    else if(i==m_y*(rank+1)-1 && i!=m)
		      {
			average = (V[up] + V[left] + V[right] + buff_2[j]);
			V[transformer(i,j)] = average;
		      }
		    else
		      {
			if (!(j >= x0 && j <= x1 && i == y0) && !(j >= x0 && j <= x1 && i == y1))
			  {
			    average = (V[up] + V[down] + V[left] + V[right])/4;
			    V[transformer(i,j)] = average;
			  }
		      }
		  }
		else if(rank==size-1)
		  {
		    if(i==m)
		      {
			average = 0;
			V[transformer(i,j)] = average;
		      }
		    else if(i==m_y*rank && i!=0)
		      {
			average = (V[down] + V[left] + V[right] + buff_1[j])/4;
			V[transformer(i,j)] = average;
		      }
		    else
		      {
			if (!(j >= x0 && j <= x1 && i == y0) && !(j >= x0 && j <= x1 && i == y1))
			  {
			    average = (V[up] + V[down] + V[left] + V[right])/4;
			    V[transformer(i,j)] = average;
			  }
		      }
		  }
		else
		  {
		    if(i==m_y*rank)
		      {
			average = (V[down] + V[left] + V[right] + buff_1[j])/4;
			V[transformer(i,j)] = average;
		      }
		    else if(i==m_y*(rank+1)-1)
		      {
			average = (V[up] + V[down] + V[left] + V[right])/4;
			V[transformer(i,j)] = average;
		      }
		    else
		      {
			if (!(j >= x0 && j <= x1 && i == y0) && !(j >= x0 && j <= x1 && i == y1))
			  {
			    average = (V[up] + V[down] + V[left] + V[right])/4;
			    V[transformer(i,j)] = average;
			  }
		      }
		  }
	    }
	}
	
        n += 1;
        
        for(k = 0; k < m; k++)
        {
            // asignar los valores a enviar de frontera
	    front_1[k] = V[transformer(m_y*rank, k)];
            front_2[k] = V[transformer(m_y*(rank +1), k)];
        }
        
        /*
         * ENVIO DE FRONTERAS
	 Creo que es posible que las fronteras se actualicen antes de ser enviadas, preferiría enviar una copia y recibir en las variables front, o enviar las variables front y recibir en una copia.
         * */
        if (rank == 0)
        {
            MPI_Irecv(buff_2, m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Isend(front_2, m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &send_request);
        }
        else if (rank == size-1)
        {
            MPI_Irecv(buff_1, m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Isend(front_1, m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &send_request);
        }
        else
        {
            MPI_Irecv(buff_1, m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Isend(front_1, m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &send_request);
            MPI_Irecv(buff_2, m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Isend(front_2, m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &send_request);
        }
        for(k = 0; k < m; k++)
        {
            //asignar las fronteras al potencial
	  //Creo que está bien
            V[transformer(m_y*rank, k)] = buff_1[k];
            V[transformer(m_y*(rank +1), k)] = buff_2[k];
        }
    }
	
    double *complete = NULL;
    if(rank == 0)
    {
        complete = (double *)malloc((size*m*m_y)*sizeof(double));
    }
      
    MPI_Gather(V, m*m_y, MPI_DOUBLE, complete, m*m_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        for(i=0; i < m; i++)
        {
            for(j=0; j < m; j++)
            {
                printf("%f\n", complete[i*m + j]);
            }
        }
    }
    MPI_Finalize();    
}

int transformer(int i, int j)
{	
	return (i - m_y*rank)*m + j;
}

double *init(int x0, int x1, int y0, int y1, double *array)
{	
	int a;
    if (y0 > min_range && y0 < max_range)
    {
        for(a = x0; a <= x1; a++)
        {
            array[transformer(y0, a)] = V0/2;
        }
    }
    else if (y1 > min_range && y1 < max_range)
    {
        for(a = x0; a <= x1; a++)
        {
            array[transformer(y1, a)] = -V0/2;
        }
    }
	return array;
}
