#include "Communicator.h"

#include <iostream>
#include <cmath>
#include <cassert>

#include <netdb.h>
#include <unistd.h>

static int send_all(int fd, void *buffer, unsigned int size);
static int recv_all(int fd, void *buffer, unsigned int size);

void Communicator::sendState(int agentId, int info, std::vector<double>& state, double reward)
{
    //std::cout <<"Send: "<<agentId<<" "<<msgID++<<" "<< info<<" ";
    for (int j=0; j<nStates; j++) {
        if (std::isnan(state[j]) || std::isinf(state[j]) ) abort();
    }

    {int *ptr=(int*)(dataout);   *ptr=agentId;}
    {int *ptr=(int*)(dataout+1); *ptr=info;  }

    for (int j=0; j<nStates; j++) {
	assert(j+2<sizeout);
        *(dataout +j+2) = state[j];
        //std::cout << state[j] << " " << " (" << j << ") ";
    }

    *(dataout +2+nStates) = reward;
    //std::cout << reward << "\n";
    //fflush(0);
    printf("Sending.\n"); fflush(0);
    send_all(sock, dataout, sizeout);
    printf("Sent.\n"); fflush(0);
    if (info == 2)  msgID = 0;
}

void Communicator::recvAction(std::vector<double>& actions)
{
    //fflush(0);
    int bytes = recv_all(sock, datain, sizein);
    if (bytes <= 0) {
        printf("selectserver: socket hung up\n");
        fflush(0);
        abort();
    }
    std::cout << "Recv: ";
    for (int j=0; j<nActions; j++) {
        actions[j] = *(datain +j);
        std::cout << actions[j] << " ";
    }
    std::cout << "\n";
}

Communicator::~Communicator()
{
    std::cout<<o.str()<<std::endl;
    close(sock);
    free(datain);
    free(dataout);
}

Communicator::Communicator(int _sockID, int _statedim, int _actdim): workerid(_sockID), nActions(_actdim), nStates(_statedim), sock(0), sizein(0), msgID(0)
{
    sprintf(SOCK_PATH, "%s%d", "/tmp/smarties_sock_", workerid);
    printf("SOCK_PATH=->%s<-\n", SOCK_PATH);

    sizeout = (3+nStates)*sizeof(double);
    dataout = (double *) malloc(sizeout);
    memset(dataout, 0, sizeout);
    sizein  =    nActions*sizeof(double);
    datain  = (double *) malloc(sizein);
    memset(datain, 0, sizein);

    /* Create a socket */
    printf("problem dim = %d %d %d \n", nStates, sizein, sizeout);

    if ((ListenerSocket = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
        perror("socket");
        exit(1);
    }
    unlink(SOCK_PATH);
    printf("Created socket over %d.\n", _sockID); fflush(0);

    bzero(&serverAddress, sizeof(serverAddress));
    serverAddress.sun_family = AF_UNIX;
    strcpy(serverAddress.sun_path, SOCK_PATH);
    const int servlen = sizeof(serverAddress.sun_family) + strlen(serverAddress.sun_path);

    if (bind(ListenerSocket, (struct sockaddr *)&serverAddress, servlen) < 0) {
        perror("bind");
        exit(1);
    }

    /* listen (only 1)*/
    if (listen(ListenerSocket, 1) == -1) {
        perror("listen");
        exit(1);
    }

    unsigned int addr_len = sizeof(clientAddress);
    if ((sock = accept(ListenerSocket, (struct sockaddr*)&clientAddress, &addr_len)) == -1) {
        perror("accept");
        return;
    }
    else printf("selectserver: new connection from on socket %d\n", sock);
    fflush(0);
}

void Communicator::dbg(double *x, int *pn)
{
    int i, n = *pn;
    int me = getpid();	/* spanwer_id : workerid */

    printf("spanwer(%d): running task with params (", me);
    for (i = 0; i < n-1; i++)
        printf("%.6lf,", x[i]);
    printf("%.6lf)\n", x[i]);
    fflush(0);
}




/*************************************************************************/
/**************************   HELPER ROUTINES   **************************/
/*************************************************************************/

static int recv_all(int fd, void *buffer, unsigned int size)
{
    int result;
    unsigned int s=size;
    char *pos = (char*)buffer;


    do {
        result=recv(fd,pos,s,0);
        if((result!=-1)&&(result>0)) {
            s -= result;
            pos += result;
        }
        else
            return result; /*-1;*/
    } while (s>0);
    //printf("recver %f\n",*((double*)buffer));
    return size;
}

static int send_all(int fd, void *buffer, unsigned int size)
{
    int result;
    unsigned int s=size;
    char *pos = (char*)buffer;

    do {
        result=send(fd,pos,s,0);
        if((result!=-1)&&(result>0)) {
            s -= result;
            pos += result;
        }
        else return result; /*-1;*/
    } while (s>0);
    return size;
}
