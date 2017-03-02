/*
 *  main.cpp
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <mpi.h>
#include <stdlib.h>


int get_rank(MPI_Comm comm)
{
	if (comm == MPI_COMM_NULL) return -1;

	int rank;
	MPI_Comm_rank(comm, &rank);
	return rank;
}

int get_nprocs(MPI_Comm comm)
{
	if (comm == MPI_COMM_NULL) return -1;

	int nprocs;
	MPI_Comm_size(comm, &nprocs);
	return nprocs;
}


extern void server_init();
extern void server_shutdown();
extern int server_fetch_and_add();
extern int jobs_init();
extern void job_run(int, MPI_Comm);

int main(int argc, char *argv[])
{
	int rank, nprocs;

	int provided;
	MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &provided);	// peh: dropped C++ bindings
	rank = get_rank(MPI_COMM_WORLD);
	nprocs = get_nprocs(MPI_COMM_WORLD);

	int ngroups = 1;
	int group_size = nprocs;

	FILE *fp = fopen("mpi.conf", "r");
	if (fp != NULL)
	{
		fscanf(fp, "%d", &ngroups);
		fscanf(fp, "%d", &group_size);
		fclose(fp);
	}

	if (nprocs != ngroups*group_size) exit(1);

	int group_ranks[group_size];

	int group_id = rank / group_size;
	int group_base = group_id * group_size;

	for (int i = 0; i < group_size; i++)
		group_ranks[i] = group_base + i;

	MPI_Group orig_group, app_group;

	/* Extract the original group handle */
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

	MPI_Group_incl(orig_group, group_size, group_ranks, &app_group);

	MPI_Comm app_comm = MPI_COMM_NULL;
	MPI_Comm_create(MPI_COMM_WORLD, app_group, &app_comm);

	int nsims = jobs_init();

	if (get_rank(MPI_COMM_WORLD) == 0)
	{
		printf("nsims = %d\n", nsims);
	}

	/// PARALLEL LOOP
	//server_init();
	for (int i = 0; i < nsims; i++)
        {
		if ((i % ngroups) != group_id) continue;
		job_run(i, app_comm);
	}

	MPI_Barrier(MPI_COMM_WORLD);	
	//server_shutdown();
	MPI_Finalize();
	return 0;

}


