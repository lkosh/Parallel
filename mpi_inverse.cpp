#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
//#include <curses.h>
#include <time.h>
#include <omp.h>
//#include <windows.h>
#include <vector>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include "mpi.h" 

using namespace std;
int main(int argc, char **argv) //метод Гауса
{
	MPI_Init(&argc, &argv);
    cout<<endl<<endl<<"Mpi realisation of Gauss method"<<endl;
    int i=0,j=0,k=0;
    int n;
    cout<<"Enter n: ";
    cin>>n;
    int temp;
    double det = 1;
 
    const double EPS = 1E-9;
    int timein, timeout, timeres = 0;
  
    int rank, nprocs;
    int map[n];
	MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* get number of processes */

    //////////////////////////////////////////////
    //n=2;
    srand(1);
    vector <vector<double> > Matrix (n, vector<double> (n));
    double c[n];
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            Matrix[i][j] = rand()%10;            
        }
    }
    
    
     for (i = 0; i < n ; i ++){
		for (j = 0; j < n; j ++)
			cout<<Matrix[i][j]<<" ";
		cout<<endl;
	}
 
    vector <vector<double> > E (n, vector<double> (n)); 
 
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            E[i][j] = 0.0;
 
            if (i == j)
                E[i][j] = 1.0;
        }
 
    double wall_timer = omp_get_wtime();
////////////////////////////////////////////////////////////////////
	

    MPI_Bcast (&Matrix[0][0],n*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    for(i=0; i<n; i++)
    {
        map[i]= i % nprocs;
    } 
    
    for(k=0;k<n;k++)
    {
		   MPI_Bcast (&Matrix[k][k], n-k, MPI_DOUBLE, map[k], MPI_COMM_WORLD);
		   MPI_Bcast (&E[k][k], n-k, MPI_DOUBLE,map[k], MPI_COMM_WORLD);
			for (i = 0; i < n; i ++){
				if (map[i] == rank){
					Matrix[k][i]   /= Matrix[k][k];
					E[k][i] /= Matrix[k][k];
				}
			}        
			// Идём по строкам, которые расположены ниже исходной
			for (int i = k + 1; i < n; ++i)
			{
				// Запоминаем множитель - элемент очередной строки,
				// расположенный под диагональным элементом исходной
				// строки
				double multi = Matrix[i][k];
				
				// Отнимаем от очередной строки исходную, умноженную
				// на сохранённый ранее множитель как в исходной,
				// так и в единичной матрице
				if (map[i] == rank){
					for (int j = 0; j < n; ++j){
						Matrix[i][j]  -= multi * Matrix[k][j];
						E[i][j] -= multi * E[k][j];
					}
				}
			}
        //for(i= k+1; i<n; i++) 
        //{
            //if(map[i] == rank)
            //{
                //c[i]=Matrix[i][k]/Matrix[k][k];
            //}
        //}               
        //for(i= k+1; i<n; i++) 
        //{       
            //if(map[i] == rank)
            //{
                //for(j=0;j<n;j++)
                
                    //Matrix[i][j]=Matrix[i][j]-( c[i]*Matrix[k][j] );
                
                
            //}
        //}
    }
    cout<<"OUTPUT"<<endl;
    
    for (i = 0; i < n ; i ++){
		for (j = 0; j < n; j ++){
			if (Matrix[i][j] < 1e-14)
				Matrix[i][j] = 0;
			cout<<Matrix[i][j]<<" ";
		}
		cout<<endl;
	}
	
	
cout<< " time on wall: " <<  omp_get_wtime() - wall_timer << "\n";
}


