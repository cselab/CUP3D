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
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
//////////////////////////////////////////////////
static int get_rank(MPI_Comm comm)
{
	if (comm == MPI_COMM_NULL) return -1;

	int rank;
	MPI_Comm_rank(comm, &rank);
	return rank;
}

static int get_nprocs(MPI_Comm comm)
{
	if (comm == MPI_COMM_NULL) return -1;

	int nprocs;
	MPI_Comm_size(comm, &nprocs);
	return nprocs;
}

////////
int fd;
fpos_t pos;

void redirect_stdout_init()
{
	fflush(stdout);
	fgetpos(stdout, &pos);
	fd = dup(fileno(stdout));
	freopen("stdout.out", "w", stdout);
}


void redirect_stdout_finalize()
{
	dup2(fd, fileno(stdout));
	close(fd);
	clearerr(stdout);
	fsetpos(stdout, &pos);        /* for C9X */
}

//////
int parse2(char *line, char **argv)
{
	int argc = 0;

	while (*line != '\0') {         /* if not the end of line ....... */
		while (*line == ' ' || *line == '\t' || *line == '\n')
			*line++ = '\0';         /* replace white spaces with 0 */
		*argv++ = line;         /* save the argument position */

		if (*line != '\0' && *line != ' ' && *line != '\t' && *line != '\n')
			argc++;

		while (*line != '\0' && *line != ' ' &&
			*line != '\t' && *line != '\n')
			line++; /* skip the argument until ...*/
	}
	*argv = '\0';   /* mark the end of argument list */

	return argc;
}

/////
int *addr = NULL;
MPI_Win win;
int done = 0;
pthread_t thread;

void *server_counter(void *arg)
{
	printf("server thread...\n"); fflush(0);
	while (!done)
	{
		int msg;
		MPI_Status status;
		MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, 998, MPI_COMM_WORLD, &status);
		if (done || (msg ==33)) break;
		MPI_Send(addr, 1, MPI_INT, status.MPI_SOURCE, 999, MPI_COMM_WORLD);
		*addr += 1;
	}
	return 0;
}

void server_init()
{
	// creation of the "shared" counter on process of rank 0
	int winSz = 0;
	if (get_rank(MPI_COMM_WORLD)==0)
	{
		MPI_Alloc_mem(winSz, MPI_INFO_NULL, &addr);
		*addr = 0; // initialised to 0 since MPI_Fetch_and_op returns value *before* increment
		pthread_create(&thread, NULL, server_counter, NULL);
	}
}

int server_fetch_and_add()
{
	// atomic incrementation of the counter
	int counter, one = 1;

	MPI_Status status;
	MPI_Send(&one, 1, MPI_INT, 0, 998, MPI_COMM_WORLD);
	MPI_Recv(&counter, 1, MPI_INT, 0, 999, MPI_COMM_WORLD, &status);

	printf("fetch_and_add returns %d\n", counter); fflush(0);
	return counter;
}

void server_shutdown()
{
	if (get_rank(MPI_COMM_WORLD)==0) {
		done = 1;
		int msg = 33;
		MPI_Send(&msg, 1, MPI_INT, 0, 998, MPI_COMM_WORLD);
//		pthread_join(thread, NULL);
	}
}



//////////////////
typedef struct job_s
{
	char str[1024];
	char dir[64];
	int nnodes, ncores, nranks, ntaskspernode, nthreads;
} job_t;

void job_print(int id, job_t *job)
{
	printf("[JOB:%d] dir=%s, nnodes=%d, ncores=%d, nranks=%d, ntaskspernode=%d, nthreads=%d\n",
		id, job->dir, job->nnodes, job->ncores, job->nranks, job->ntaskspernode, job->nthreads); fflush(0);
}

#define EXECNAME	"mpcf-cluster"
#define CMDFILENAME	"queue.dat"
#define MAXJOBS	64
job_t jobs[MAXJOBS];

int jobs_init()
{
	FILE *cmdfp;
	cmdfp = fopen(CMDFILENAME, "r");
	char line[1024];

	int njobs;

	char execname[256];
	strcpy(execname, EXECNAME);

	njobs = 0;
	while (fgets(line, 1024, cmdfp)!=NULL)
	{
		if (line[0]=='#') continue;
		if (strstr(line, execname) != NULL) {
			strcpy(jobs[njobs].str, line);
			sscanf(line, "%s %d %d %d %d %d %*s", &jobs[njobs].dir, &jobs[njobs].nnodes,
					&jobs[njobs].ncores, &jobs[njobs].nranks, &jobs[njobs].ntaskspernode, &jobs[njobs].nthreads);
			njobs++;
		}
	}
	fclose(cmdfp);

	return njobs;
}


extern void cubism_main(MPI_Comm comm, int argc, char **argv);

void job_run(int i, MPI_Comm app_comm)
{
	// build argc, argv from line
	int largc;
	char *largv[256];
	char *cmdline;

	cmdline = strstr(jobs[i].str, EXECNAME);
	largc = parse2(cmdline, largv);

	if (get_rank(app_comm) == 0)
	{
		sleep(3);
		job_print(i, &jobs[i]);
//		for (int j = 0; j < largc; j++)
//			printf("arg %d: %s\n", j, largv[j]);
	}

#if 1
	char initd[256];
	char newd[256];

	getcwd(initd,256);
	sprintf(newd,"%s/%s", initd, jobs[i].dir);
	chdir(newd);	// go to the task private directory

	redirect_stdout_init();
	cubism_main(app_comm, largc, largv);
	redirect_stdout_finalize();

	chdir(initd);	// go up one level
#endif

}
