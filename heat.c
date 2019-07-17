#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define m 20 
#define n 20 
#define epsilon    0.001   
#define tolerance  0.001
#define max(a,b)   a > b ? a : b

int main(int argc, char** argv) {
	FILE* fp = fopen("dung.txt", "w");
	int i, j, rank, NP, mc;
	float delta_local, delta_global, temp;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NP);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	mc = n / NP;
	float *C = (float *)malloc(sizeof(float)*m*n);
	float *Cs = (float *)malloc(mc*n * sizeof(float));
	float *Cu = (float *)malloc(n * sizeof(float));
	float *Cd = (float *)malloc(n * sizeof(float));

	if (rank == 0) {
		KhoiTao(C);

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (j < n - 1) {
					fprintf(fp, "%.2f,", *(C + i * n + j));
				}
				else {
					fprintf(fp, "%.2f\n", *(C + i * n + j));
				}
			}
		}
	}

	//Devide Input And Send To All CPU
	MPI_Scatter(C, mc*n, MPI_FLOAT, Cs, mc*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	//Red-Black Ordering
	int iterator = 0;
	do {
		delta_global = 0.00;
		delta_local = 0.00;

		//Exchange boundary strips with neighboring processors;
		//Exchange Cu
		if (rank == 0) {
			MPI_Send(Cs + (mc - 1)*n, n, MPI_FLOAT, rank + 1, rank, MPI_COMM_WORLD);
		}
		else if (rank == NP - 1) {
			MPI_Recv(Cu, n, MPI_FLOAT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
		}
		else {
			MPI_Send(Cs + (mc - 1)*n, n, MPI_FLOAT, rank + 1, rank, MPI_COMM_WORLD);
			MPI_Recv(Cu, n, MPI_FLOAT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
		}
		//Exchange Cd
		if (rank == NP - 1) {
			MPI_Send(Cs, n, MPI_FLOAT, rank - 1, rank, MPI_COMM_WORLD);
		}
		else if (rank == 0) {
			MPI_Recv(Cd, n, MPI_FLOAT, rank + 1, rank + 1, MPI_COMM_WORLD, &status);
		}
		else {
			MPI_Send(Cs, n, MPI_FLOAT, rank - 1, rank, MPI_COMM_WORLD);
			MPI_Recv(Cd, n, MPI_FLOAT, rank + 1, rank + 1, MPI_COMM_WORLD, &status);
		}

		//Update Red Points In This Processor
		for (i = 0; i < mc; i++)
			for (j = 0; j < n; j++)
				if ((i + j) % 2 == 0) {
					//Red Points
					temp = *(Cs + n * i + j);
					if ((rank == 0 && i == 0) || (rank == NP - 1 && i == mc - 1))
						//Global Boundary Condition
						*(Cs + n * i + j) = 25;
					else {
						//Update Red Point Cs(i,j)
						//Local Boundary 
						float up    = i == 0      ? *(Cu + j)             : *(Cs + n * (i - 1) + j);
						float down  = i == mc - 1 ? *(Cd + j)             : *(Cs + n * (i + 1) + j);
						//Circulation
						float left  = j == 0      ? *(Cs + n * i + n - 1) : *(Cs + n * i + (j - 1));
						float right = j == n - 1  ? *(Cs + n * i + 0)     : *(Cs + n * i + (j + 1));
						*(Cs + n * i + j) = 0.25*(up + down + left + right);
					}
					//Calculate Local Delta 
					delta_local = max(delta_local, fabs(*(Cs + n * i + j) - temp));
				}

		//Exchange boundary strips with neighboring processors;
		//Exchange Top
		if (rank == 0) {
			MPI_Send(Cs + (mc - 1)*n, n, MPI_FLOAT, rank + 1, rank, MPI_COMM_WORLD);
		}
		else if (rank == NP - 1) {
			MPI_Recv(Cu, n, MPI_FLOAT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
		}
		else {
			MPI_Send(Cs + (mc - 1)*n, n, MPI_FLOAT, rank + 1, rank, MPI_COMM_WORLD);
			MPI_Recv(Cu, n, MPI_FLOAT, rank - 1, rank - 1, MPI_COMM_WORLD, &status);
		}
		//Exchange Bot
		if (rank == NP - 1) {
			MPI_Send(Cs, n, MPI_FLOAT, rank - 1, rank, MPI_COMM_WORLD);
		}
		else if (rank == 0) {
			MPI_Recv(Cd, n, MPI_FLOAT, rank + 1, rank + 1, MPI_COMM_WORLD, &status);
		}
		else {
			MPI_Send(Cs, n, MPI_FLOAT, rank - 1, rank, MPI_COMM_WORLD);
			MPI_Recv(Cd, n, MPI_FLOAT, rank + 1, rank + 1, MPI_COMM_WORLD, &status);
		}

		//Update Black Points In This Processor
		for (i = 0; i < mc; i++)
			for (j = 0; j < n; j++)
				if ((i + j) % 2 != 0) {
					//Black Points
					temp = *(Cs + n * i + j);
					if ((rank == 0 && i == 0) || (rank == NP - 1 && i == mc - 1))
						//Global Boundary Condition
						*(Cs + n * i + j) = 25;
					else {
						//Update Black Point
						//Local Boundary
						float up    = i == 0      ? *(Cu + j)             : *(Cs + n * (i - 1) + j);
						float down  = i == mc - 1 ? *(Cd + j)             : *(Cs + n * (i + 1) + j);
						//Circulation
						float left  = j == 0      ? *(Cs + n * i + n - 1) : *(Cs + n * i + (j - 1));
						float right = j == n - 1  ? *(Cs + n * i + 0)     : *(Cs + n * i + (j + 1));
						*(Cs + n * i + j) = 0.25*(up + down + left + right);
					}
					//Calculate Delta
					delta_local = max(delta_local, fabs(*(Cs + n * i + j) - temp));
				}

		MPI_Allreduce(&delta_local, &delta_global, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

		//Send Output To Root
		MPI_Gather(Cs, mc*n, MPI_FLOAT, C, mc*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("Iterator: %d, delta: %.6f \n", iterator, delta_global);
			iterator++;

			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++)
					printf("%.2f ", *(C + n * i + j));
				printf("\n");
			}

			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					if (j < n - 1) {
						fprintf(fp, "%.2f,", *(C + i * n + j));
					}
					else {
						fprintf(fp, "%.2f\n", *(C + i * n + j));
					}
				}
			}

		}



	} while (delta_global > tolerance);

	MPI_Finalize();
	return 0;
}

void KhoiTao(float *C) {
	int i, j;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (i >= (m / 2 - 5) && i < (m / 2 + 5) 
			 && j >= (n / 2 - 5) && j < (m / 2 + 5))
				*(C + n * i + j) = 80;
			else
				*(C + n * i + j) = 25;
}




