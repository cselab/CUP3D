#include "communicator.h"
#include <iostream>
#include <cmath>
#include "Util/util.h"

void Communicator::sendState(int agentId, int info, std::vector<double>& state, double reward)
{
    o <<"Send: "<<agentId<<" "<<msgID++<<" "<< info<<" ";
    for (int j=0; j<nStates; j++) {
        if (std::isnan(state[j]) || std::isinf(state[j]) ) abort();
    }
    
    {int *ptr=(int*)(dataout);   *ptr=agentId;}
    {int *ptr=(int*)(dataout+1); *ptr=info;  }

    for (int j=0; j<nStates; j++) {
        *(dataout +j+2) = state[j];
        o << state[j] << " ";
    }
    
    *(dataout +2+nStates) = reward;
    o << reward << "\n";
    
    send_all(sock, dataout, sizeout);
    if (info == 2)  msgID = 0;
}
    
void Communicator::recvAction(std::vector<double>& actions)
{
    int bytes = recv_all(sock, datain, sizein);
    if (bytes <= 0) {
        printf("selectserver: socket hung up\n");
        fflush(0);
        abort();
    }
    o << "Recv: ";
    for (int j=0; j<nActions; j++) {
        actions[j] = *(datain +j);
        o << actions[j] << " ";
    }
    o << "\n";
}
    
Communicator::~Communicator()
{
    std::cout<<o.str()<<std::endl;
    close(sock);
}
    
Communicator::Communicator(int _sockID, int _statedim, int _actdim): workerid(_sockID), nActions(_actdim), nStates(_statedim), sock(0), sizein(0), msgID(0)
{
    sprintf(SOCK_PATH, "%s%d", "/tmp/sock_", workerid);
    printf("SOCK_PATH=->%s<-\n", SOCK_PATH);

    sizeout = (3+nStates)*sizeof(double);
    dataout = (double *) malloc(sizeout);
    sizein  =    nActions*sizeof(double);
    datain  = (double *) malloc(sizein);
    
    /* Create a socket */
    printf("problem dim = %i %d %d \n", nStates, sizein, sizeout);
    
    if ((ListenerSocket = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
        perror("socket");
        exit(1);
    }
    unlink(SOCK_PATH);
    
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
