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
//    printf("%d %d\n", size, rank);
    
    //MPI_Request send_request, recv_request;
    //MPI_Status status;
    
	int up, down, left, right, j=1, k = 0, n = 0;
    int x0, x1, y0, y1, N, i, until;
    double average, der, h = (double)L/m;
	
    m_y = m/size;
    N = 2*m*m_y;
    max_range = m_y*(rank + 1);
    min_range = m_y*(rank);
        
	x0 = m/2 - l/(h*2) - 1;
	x1 = m/2 + l/(h*2) - 1;        
	y0 = m/2 - d/(h*2) - 1;
	y1 = m/2 + d/(h*2) - 1;
    
    
    
    double *V = malloc(m*(m_y+1)*sizeof(double));
    double *buffer_1 = malloc(m*sizeof(double));
    double *buffer_2 = malloc(m*sizeof(double));
    double *front_1 = malloc(m*sizeof(double));
    double *front_2 = malloc(m*sizeof(double));
    V = init(x0, x1, y0, y1, V);

    while (n < N)
    {	
        i = m_y*rank + 1;
        until = (rank + 1)*m_y;
        for(i ; i < until; i++)
        {
            for(j=1;j < m-1; j++)
            {
                up = transformer(i-1, j);
                down = transformer(i+1, j);
                left = transformer(i, j-1);
                right = transformer(i, j+1);
                if (!(j >= x0 && j <= x1 && i == y0) && !(j >= x0 && j <= x1 && i == y1))
                {
                    average = (V[up] + V[down] + V[left] + V[right])/4;
                    V[transformer(i,j)] = average;
                }
            }
        }
        n += 1;
        
        if(rank == 0)
        {        
            for(k = 0; k < m; k++)
            {
                // asignar los valores a enviar de frontera
                front_2[k] = V[transformer(until-1, k)];
            }
        }
        if(rank == size - 1)
        {        
            for(k = 0; k < m; k++)
            {
                // asignar los valores a enviar de frontera
                front_1[k] = V[transformer(m_y*rank+1, k)];
            }
        }
        else
        {
            for(k = 0; k < m; k++)
            {
                front_1[k] = V[transformer(m_y*rank+1, k)];
                front_2[k] = V[transformer(until-1, k)];
            }
        }
        
        /*
         * ENVIO DE FRONTERAS 
         * */
        if (rank == 0)
        {
            MPI_Send(front_2, m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(buffer_2, m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if(rank == size - 1)
        {
            MPI_Recv(buffer_1, m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(front_1, m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Recv(buffer_1, m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(front_1, m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Send(front_2, m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(buffer_2, m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        }
        for(k = 0; k < m; k++)
        {
            //asignar las fronteras al potencial
            V[transformer(m_y*rank, k)] = buffer_1[k];
            V[transformer(until, k)] = buffer_2[k];
        }
    }
	
    double *complete = NULL;
    if(rank == 0)
    {
        complete = (double *)malloc((size*m*m_y)*sizeof(double));
    }
      
    MPI_Gather(V, m*m_y, MPI_DOUBLE, complete, m*m_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(rank==0)
      {
    int sizes = 256;
    FILE *f = fopen("potential.dat","w");
  for(i=0;i<sizes;i++)
    {
      for(j=0;j<sizes;j++)
	{
	  fprintf(f,"%f ",complete[sizes*i+j]);
	}
      fprintf(f,"\n");
    }
  fclose(f);

  FILE *fieldx = fopen("fieldx.dat","w");
  for(i=0;i<sizes;i++)
    {
      for(j=0;j<sizes;j++)
	{
	  if(j==0)
	    {
	      der = -(complete[sizes*i+j+1]-complete[sizes*i+j])/h;
	      fprintf(fieldx,"%f ",der);
	    }
	  else if(j==size-1)
	    {
	      der = -(complete[sizes*i+j]-complete[sizes*i+j-1])/h;
	      fprintf(fieldx,"%f ",der);
	    }
	  else
	    {
	      der = -(complete[sizes*i+j+1]-complete[sizes*i+j-1])/(2.0*h);
	      fprintf(fieldx,"%f ",der);
	    }
	}
      fprintf(fieldx,"\n");
    }
  fclose(fieldx);

  FILE *fieldy = fopen("fieldy.dat","w");
  for(i=0;i<sizes;i++)
    {
      for(j=0;j<sizes;j++)
	{
	  if(i==0)
	    {
	      der = (complete[sizes*(i+1)+j]-complete[sizes*i+j])/h;
	      fprintf(fieldy,"%f ",der);
	    }
	  else if(i==size-1)
	    {
	      der = (complete[sizes*i+j]-complete[sizes*(i-1)+j])/h;
	      fprintf(fieldy,"%f ",der);
	    }
	  else
	    {
	      der = (complete[sizes*(i+1)+j]-complete[sizes*(i-1)+j])/(2.0*h);
	      fprintf(fieldy,"%f ",der);
	    }
	}
      fprintf(fieldy,"\n");
    }
    fclose(fieldy);
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
