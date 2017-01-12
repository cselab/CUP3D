#pragma once
#include <vector>
#include <sstream>
#include <sys/un.h>

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
