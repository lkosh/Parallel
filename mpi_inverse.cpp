#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
//#include <curses.h>
#include <time.h>
#include <omp.h>
//#include <windows.h>

#include <math.h>
#include <iomanip>
#include <cstdlib>
#include "mpi.h" 

using namespace std;
void PrintM (double **Matrix , int n ){
	 for (int i = 0; i < n ; i ++){
			 cout<<"->  ";
			for (int j = 0; j < n; j ++)
				cout<<Matrix[i][j]<<" ";
			cout<<endl;
		}
}
int main(int argc, char **argv) //метод Гауса
{
	
    //cout<<endl<<endl<<"Mpi realisation of Gauss method"<<endl;
    int i=0,j=0,k=0;
    int n = 5;
    //cout<<"Enter n: ";
    //cin>>n;
    int temp;
    double det = 1;
 
    const double EPS = 1E-9;
    int timein, timeout, timeres = 0;
  
    int rank=0, nprocs;
    int map[n];
   
	
	MPI_Init(&argc, &argv);
	MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* get current process id */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* get number of processes */	
    //////////////////////////////////////////////
    //n=2;
    srand(1);
    //cout<<"M start"<<endl;
	double Matrix[n*n];
	double E[n*n];
    double c[n];
	if (!rank){
		for (int i=0;i<n;i++)
		{
			for (int j=0;j<n;j++)
			{
				Matrix[i*n +j] = rand()%10;  
				//cout<<i<<" "<<j<<endl;          
			}
		}
	//cout<<"M "<<endl;
	   
		
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				E[i*n + j] = 0.0;
	 
				if (i == j)
					E[i*n + j] = 1.0;
			}
	}
	//cout<<"M end"<<endl;
    //if (!rank)
		//cout<<nprocs<<endl;
	//cout<<rank;
	//if (!rank){
		 //for (i = 0; i < n ; i ++){
			 //cout<<"->  ";
			//for (j = 0; j < n; j ++)
				//cout<<Matrix[i][j]<<" ";
			//cout<<endl;
		//}
	//}
	//if (!rank)
	double wall_timer = MPI_Wtime();//omp_get_wtime();
////////////////////////////////////////////////////////////////////
	
	//cout<<"start"<<endl;
    MPI_Bcast (&Matrix[0],n*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast (&E[0],n*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(i=0; i<n; i++)
    {
        map[i]= i % nprocs;
		
    } 
  //  cout<<"start 2 "<<endl;
    double tmp;
    for(k=0;k<n;k++)
    {
		 MPI_Bcast (&Matrix[k*n + k], n-k, MPI_DOUBLE, map[k], MPI_COMM_WORLD);
		 MPI_Bcast (&E[k * n + k], n-k, MPI_DOUBLE,map[k], MPI_COMM_WORLD);
		if (fabs(Matrix[k*n +k]) < 1e-8 && !rank)
        {
            // Ключ, говорязий о том, что был произведён обмен строк
            bool changed = false;
            
            // Идём по строкам, расположенным ниже исходной
            for (int i = k + 1; i < n; ++i)
            {
                // Если нашли строку, где в том же столбце
                // имеется ненулевой элемент
                if (fabs(Matrix[i*n + k]) > 1e-8)
                {
                    // Меняем найденную и исходную строки местами
                    // как в исходной матрице, так и в единичной
                    for (int h = 0; h < n; h ++){
						tmp = Matrix[k*n + h];
						Matrix[k*n + h] = Matrix[i*n + h];
						Matrix[i*n + h] = tmp;
						
						tmp = E[k*n + h];
						E[k*n + h] = E[i*n + h];
						E[i*n + h] = tmp;
                    }
                    // Взводим ключ - сообщаем о произведённом обмене строк
                    changed = true;
                    
                    break;
                }
            }
          //  cout<<"sth"<<endl;
            // Если обмен строк произведён не был - матрица не может быть
            // обращена
            if (!changed)
            {
                cout<<"Matrix cannot inverse!!!"<<endl;
                // Сообщаем о неудаче обращения
            }
        }  
		  
			for (j = 0; j < n; j ++){
				if (map[j] == rank){
					
					Matrix[k*n + j]   /= Matrix[k*n + k];
					E[k*n + j] /= Matrix[k*n + k];
				}
			}        
		
			
			// Идём по строкам, которые расположены ниже исходной
			for (int i = k + 1; i < n; ++i)
			{
				// Запоминаем множитель - элемент очередной строки,
				// расположенный под диагональным элементом исходной
				// строки
				double multi = Matrix[i*n + k];
				
				// Отнимаем от очередной строки исходную, умноженную
				// на сохранённый ранее множитель как в исходной,
				// так и в единичной матрице
				if (map[i] == rank){
					for (int j = 0; j < n; ++j){
						Matrix[i*n + j]  -= multi * Matrix[k*n + j];
						E[i*n + j] -= multi * E[k*n + j];
					}
				}
			}

    }
	
	 // Проходим по вернхней треугольной матрице, полученной
    // на прямом ходе, снизу вверх
    // На данном этапе происходит обратный ход, и из исходной
    // матрицы окончательно формируется единичная, а из единичной -
    // обратная
    for (int k = n - 1; k > 0; --k)
    {
        // Идём по строкам, которые расположены выше исходной
        for (int i = k - 1; i + 1 > 0; --i)
        {
            // Запоминаем множитель - элемент очередной строки,
            // расположенный над диагональным элементом исходной
            // строки
            double multi = Matrix[i*n + k];
            
            // Отнимаем от очередной строки исходную, умноженную
            // на сохранённый ранее множитель как в исходной,
            // так и в единичной матрице
            if (map[i] == rank){
				for (int j = 0; j < n; ++j)
				{
					Matrix[i*n + j]   -= multi * Matrix[k*n + j];
					E[i*n + j] -= multi * E[k*n + j];
				}
			}
        }
    }
     MPI_Barrier(MPI_COMM_WORLD);
 
    //if (!rank){
		//cout<<"OUTPUT 2"<<endl;
		//for (i = 0; i < n*n ; i ++){
			//cout<<Matrix[i]<<" ";
			//cout<<endl;
		//}
	//}
	if (rank ==0){
		cout<< " time on wall: " <<  MPI_Wtime() - wall_timer << "\n";
		cout<<" n "<< n<< endl;
	}
}


