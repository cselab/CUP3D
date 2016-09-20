#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <sys/time.h>
#include <ftw.h>
#include <sys/un.h>
#include <netdb.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
//#include <omp.h>
#include <vector>
#include <fstream>
#include <sstream>

class Communicator
{
    const int workerid, nActions, nStates;
    int sock, ListenerSocket, sizein, sizeout;
    int msgID;
    std::ostringstream o;
    char SOCK_PATH[256];
    struct sockaddr_un server_addr;
    struct sockaddr_un serverAddress;
    struct sockaddr_un clientAddress;
public:
    double *datain, *dataout;
    
    void sendState(int agentId, int info, std::vector<double> & state, double reward);
    
    void recvAction(std::vector<double> & actions);
    
    ~Communicator();
    
    Communicator(int _sockID, int _statedim, int _actdim);
    
    void dbg(double *x, int *pn);
};
