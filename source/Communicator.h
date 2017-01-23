#pragma once

#include <sstream>
#include <sys/un.h>
#include <mpi.h>
#include <vector>

#define _AGENT_STATUS int
#define _AGENT_FIRSTCOMM 1
#define _AGENT_NORMCOMM  0
#define _AGENT_LASTCOMM  2
#define _AGENT_FAILCOMM  -1

class Communicator
{
    const int workerid, nActions, nStates, isServer;
    int Socket, ListenerSocket, sizein, sizeout;
    std::ostringstream o;
    int msgID, rank_MPI, size_MPI;

    char SOCK_PATH[256];
    struct sockaddr_un serverAddress;
    struct sockaddr_un clientAddress;
    double *datain, *dataout;

public:
    MPI_Comm comm_MPI;
    void setupServer();
    void setupClient(const int iter, const std::string execpath);
    void closeSocket();

    void dbg(double *x, int *pn);

    void sendState(int agentId,
                   int info,
                   std::vector<double>& state,
                   double reward);
    void recvState(int& agentId,
                   int& info,
                   std::vector<double>& state,
                   double& reward);

    // assumption is that Action always follows a State:
    // no need to specify agent
    void recvAction(std::vector<double> & actions);
    void sendAction(std::vector<double> & actions);

    ~Communicator();
    Communicator(int _sockID, int _statedim, int _actdim, bool _isserver=0, bool _issim=1);
    Communicator(int _sockID, int _statedim, int _actdim, MPI_Comm comm);
};
