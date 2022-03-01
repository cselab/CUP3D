//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#include "Fish.h"
#include "FishLibrary.h"

#include <Cubism/ArgumentParser.h>
#include <Cubism/HDF5Dumper_MPI.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

Fish::Fish(SimulationData&s, ArgumentParser&p) : Obstacle(s, p)
{
  p.unset_strict_mode();
  Tperiod = p("-T").asDouble(1.0);
  phaseShift = p("-phi").asDouble(0.0);

  //PID knobs
  bCorrectTrajectory = p("-Correct").asBool(false);
  bCorrectPosition = p("-bCorrectPosition").asBool(false);

  // MPI datatypes (used for load-balancing when creating the fish surface)
  int array_of_blocklengths[2]       = {4, 1};
  MPI_Aint array_of_displacements[2] = {0, 4 * sizeof(double)};
  MPI_Datatype array_of_types[2]     = {MPI_DOUBLE, MPI_LONG};
  MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_BLOCKID);
  MPI_Type_commit(&MPI_BLOCKID);

  const int Z = FluidBlock::sizeZ;
  const int Y = FluidBlock::sizeY;
  const int X = FluidBlock::sizeX;
  int array_of_blocklengths1[2]       = {Z*Y*X*3 + (Z+2)*(Y+2)*(X+2), Z*Y*X};
  MPI_Aint array_of_displacements1[2] = {0, (Z*Y*X*3 + (Z+2)*(Y+2)*(X+2)) * sizeof(double)};
  MPI_Datatype array_of_types1[2]     = {MPI_DOUBLE, MPI_INT};
  MPI_Type_create_struct(2, array_of_blocklengths1, array_of_displacements1, array_of_types1, &MPI_OBSTACLE);
  MPI_Type_commit(&MPI_OBSTACLE);
}

Fish::~Fish()
{
  if(myFish not_eq nullptr) delete myFish;
  MPI_Type_free(&MPI_BLOCKID);
  MPI_Type_free(&MPI_OBSTACLE);
}

void Fish::integrateMidline()
{
  myFish->integrateLinearMomentum();

  myFish->integrateAngularMomentum(sim.dt);
}

std::vector<VolumeSegment_OBB> Fish::prepare_vSegments()
{
  /*
    - VolumeSegment_OBB's volume cannot be zero
    - therefore no VolumeSegment_OBB can be only occupied by extension midline
      points (which have width and height = 0)
    - performance of create seems to decrease if VolumeSegment_OBB are bigger
    - this is the smallest number of VolumeSegment_OBB (Nsegments) and points in
      the midline (Nm) to ensure at least one non ext. point inside all segments
   */
  const int Nsegments = std::ceil((myFish->Nm-1.)/8);
  const int Nm = myFish->Nm;
  assert((Nm-1)%Nsegments==0);

  std::vector<VolumeSegment_OBB> vSegments(Nsegments);
  #pragma omp parallel for schedule(static)
  for(int i=0; i<Nsegments; ++i)
  {
    const int nextidx = (i+1)*(Nm-1)/Nsegments;
    const int idx = i * (Nm-1)/Nsegments;
    // find bounding box based on this
    Real bbox[3][2] = {{1e9, -1e9}, {1e9, -1e9}, {1e9, -1e9}};
    for(int ss=idx; ss<=nextidx; ++ss)
    {
      const Real xBnd[4] = {myFish->rX[ss] + myFish->norX[ss]*myFish->width[ss],
                            myFish->rX[ss] - myFish->norX[ss]*myFish->width[ss],
                            myFish->rX[ss] + myFish->binX[ss]*myFish->height[ss],
                            myFish->rX[ss] - myFish->binX[ss]*myFish->height[ss]};
      const Real yBnd[4] = {myFish->rY[ss] + myFish->norY[ss]*myFish->width[ss],
                            myFish->rY[ss] - myFish->norY[ss]*myFish->width[ss],
                            myFish->rY[ss] + myFish->binY[ss]*myFish->height[ss],
                            myFish->rY[ss] - myFish->binY[ss]*myFish->height[ss]};
      const Real zBnd[4] = {myFish->rZ[ss] + myFish->norZ[ss]*myFish->width[ss],
                            myFish->rZ[ss] - myFish->norZ[ss]*myFish->width[ss],
                            myFish->rZ[ss] + myFish->binZ[ss]*myFish->height[ss],
                            myFish->rZ[ss] - myFish->binZ[ss]*myFish->height[ss]};
      const Real maxX=std::max({xBnd[0],xBnd[1],xBnd[2],xBnd[3]});
      const Real maxY=std::max({yBnd[0],yBnd[1],yBnd[2],yBnd[3]});
      const Real maxZ=std::max({zBnd[0],zBnd[1],zBnd[2],zBnd[3]});
      const Real minX=std::min({xBnd[0],xBnd[1],xBnd[2],xBnd[3]});
      const Real minY=std::min({yBnd[0],yBnd[1],yBnd[2],yBnd[3]});
      const Real minZ=std::min({zBnd[0],zBnd[1],zBnd[2],zBnd[3]});
      bbox[0][0] = std::min(bbox[0][0], minX);
      bbox[0][1] = std::max(bbox[0][1], maxX);
      bbox[1][0] = std::min(bbox[1][0], minY);
      bbox[1][1] = std::max(bbox[1][1], maxY);
      bbox[2][0] = std::min(bbox[2][0], minZ);
      bbox[2][1] = std::max(bbox[2][1], maxZ);
    }

    vSegments[i].prepare(std::make_pair(idx,nextidx), bbox, sim.hmin);
    vSegments[i].changeToComputationalFrame(position,quaternion);
  }
  return vSegments;
}

using intersect_t = std::vector<std::vector<VolumeSegment_OBB*>>;
intersect_t Fish::prepare_segPerBlock(vecsegm_t& vSegments)
{
  MyBlockIDs.clear();
  for (size_t j = 0 ; j < MySegments.size(); j++)
      MySegments[j].clear();
  MySegments.clear();

  const std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
  std::vector<std::vector<VolumeSegment_OBB*>> ret(vInfo.size());

  // clear deformation velocities
  for(auto & entry : obstacleBlocks) {
    if(entry == nullptr) continue;
    delete entry;
    entry = nullptr;
  }
  obstacleBlocks.resize(vInfo.size(), nullptr);

  //#pragma omp parallel for schedule(dynamic, 1)
  for(size_t i=0; i<vInfo.size(); ++i)
  {
    const BlockInfo & info = vInfo[i];
    const FluidBlock & b = *(FluidBlock*)info.ptrBlock;

    bool hasSegments = false;
    for(size_t s=0; s<vSegments.size(); ++s)
      if(vSegments[s].isIntersectingWithAABB(b.min_pos.data(), b.max_pos.data())) {
        VolumeSegment_OBB*const ptr  = & vSegments[s];
        ret[info.blockID].push_back( ptr );

        //#pragma omp critical
        {
          if (!hasSegments)
          {
              hasSegments = true;
              MyBlockIDs.push_back({info.h,info.origin[0],info.origin[1],info.origin[2],info.blockID});
              MySegments.resize(MySegments.size()+1);
          }
          MySegments.back().push_back(s);        
        }
      }

    // allocate new blocks if necessary
    if( ret[info.blockID].size() > 0 ) {
      assert( obstacleBlocks[info.blockID] == nullptr );
      ObstacleBlock * const block = new ObstacleBlock();
      assert(block not_eq nullptr);
      obstacleBlocks[info.blockID] = block;
      block->clear();
    }
  }
  return ret;
}

void Fish::writeSDFOnBlocks(std::vector<VolumeSegment_OBB> & vSegments)
{
#if 1 //no load-balancing here
  #pragma omp parallel
  {
    PutFishOnBlocks putfish(myFish, position, quaternion);
    #pragma omp for
    for (size_t j=0 ; j < MyBlockIDs.size(); j++)
    {
      std::vector<VolumeSegment_OBB*> S;
      for (size_t k = 0 ; k < MySegments[j].size() ; k++)
        S.push_back(& vSegments[MySegments[j][k]]);
      ObstacleBlock*const block = obstacleBlocks[MyBlockIDs[j].blockID];
      putfish(MyBlockIDs[j].h, MyBlockIDs[j].origin_x, MyBlockIDs[j].origin_y, MyBlockIDs[j].origin_z, block, S);
    }
  }
#else //load-balancing - there's a bug with message tags here!!!

  MPI_Comm comm = grid->getCartComm();
  const int rank = grid->rank();
  const int size = grid->get_world_size();
  std::vector<std::vector<int>> OtherSegments;

  //Each MPI rank owns two arrays: 
  // MyBlockIDs[]: a list of blocks that have at least one segment
  // MySegments[i][j] : a list of integers j=0,... for MyBlockIDs[i]

  //Determine the total load and how much load corresponds to each rank:
  //All N blocks are indexed for 0 to N-1. Each rank computes the range of 
  //indices that it must have, in order to have an equal load distribution.
  int b = (int) MyBlockIDs.size();
  std::vector<int> all_b(size);
  MPI_Allgather(&b, 1, MPI_INT, all_b.data(), 1, MPI_INT, comm);

  int total_load = 0;
  for (int r = 0 ; r < size ; r++) total_load+= all_b[r];
  int my_load =  total_load / size;
  if (rank < (total_load % size) ) my_load += 1;

  std::vector<int> index_start(size);
  index_start[0] = 0;
  for (int r = 1 ; r < size ; r++) index_start[r] = index_start[r-1] + all_b[r-1];

  int ideal_index = ( total_load / size ) * rank;
  ideal_index += (rank < (total_load % size)) ? rank : (total_load % size);

  //Now that each rank knows what range it should have, it will determine where to send/receive
  std::vector< std::vector<BlockID> > send_blocks(size);
  std::vector< std::vector<BlockID> > recv_blocks(size);
  for (int r = 0 ; r < size ; r ++) if (rank != r)
  {
     {  //check if I need to receive blocks
        const int a1 = ideal_index;
        const int a2 = ideal_index + my_load -1;
        const int b1 = index_start[r];
        const int b2 = index_start[r]+all_b[r]-1;
        const int c1 = max(a1,b1);
        const int c2 = min(a2,b2);
        if (c2-c1 + 1>0) recv_blocks[r].resize(c2-c1+1);
     }
     {  //check if I need to send blocks
        int other_ideal_index = ( total_load / size ) * r;
        other_ideal_index += (r < (total_load % size)) ? r : (total_load % size); 
        int other_load =  total_load / size;
        if (r < (total_load%size)) other_load += 1;
        const int a1 = other_ideal_index;
        const int a2 = other_ideal_index + other_load -1;
        const int b1 = index_start[rank];
        const int b2 = index_start[rank]+all_b[rank]-1;
        const int c1 = max(a1,b1);
        const int c2 = min(a2,b2);
        if (c2-c1 + 1>0) send_blocks[r].resize(c2-c1+1);
     }
  }

  //Send and receive the information needed to create the obstacle blocks. 
  std::vector<MPI_Request> recv_request;
  for (int r = 0 ; r < size ; r ++) if (recv_blocks[r].size() != 0)
  {
    MPI_Request req;
    recv_request.push_back(req);
    MPI_Irecv(recv_blocks[r].data(), recv_blocks[r].size(), MPI_BLOCKID, r, r*size+rank, comm, &recv_request.back());
  }
  std::vector<MPI_Request> send_request;
  int counter = 0;
  for (int r = 0 ; r < size ; r ++) if (send_blocks[r].size() != 0)
  {
    for (size_t i = 0 ; i < send_blocks[r].size() ; i ++)
    {
      send_blocks[r][i].h         = MyBlockIDs[counter + i].h        ;
      send_blocks[r][i].origin_x  = MyBlockIDs[counter + i].origin_x ;
      send_blocks[r][i].origin_y  = MyBlockIDs[counter + i].origin_y ;
      send_blocks[r][i].origin_z  = MyBlockIDs[counter + i].origin_z ;
      send_blocks[r][i].blockID   = MyBlockIDs[counter + i].blockID  ;
    }
    counter += send_blocks[r].size();
    MPI_Request req;
    send_request.push_back(req);
    MPI_Isend(send_blocks[r].data(), send_blocks[r].size(), MPI_BLOCKID, r, r +rank*size, comm, &send_request.back());
  }


  //allocate buffers for the actual data that will be sent/received
  const int sizeZ = FluidBlock::sizeZ;
  const int sizeY = FluidBlock::sizeY;
  const int sizeX = FluidBlock::sizeX;
  std::vector< std::vector<MPI_Obstacle> > send_obstacles(size);
  std::vector< std::vector<MPI_Obstacle> > recv_obstacles(size);
  for (int r = 0 ; r < size ; r++)
  {
    send_obstacles[r].resize(send_blocks[r].size());
    recv_obstacles[r].resize(recv_blocks[r].size());
  }
  MPI_Waitall(send_request.size(), send_request.data() , MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_request.size(), recv_request.data() , MPI_STATUSES_IGNORE);


  // Do the blocks I received
  for (int r = 0 ; r < size ; r++) if (recv_blocks[r].size() != 0)
  {
    for (size_t j = 0 ; j < OtherSegments.size(); j++) OtherSegments[j].clear();
    OtherSegments.clear();
    for(size_t i=0; i<recv_blocks[r].size(); ++i)
    {
      const auto & info = recv_blocks[r][i];
      bool hasSegments = false;
      for(size_t s=0; s<vSegments.size(); ++s)
      {
        double min_pos [3] = {info.origin_x + 0.5*info.h,info.origin_y + 0.5*info.h,info.origin_z + 0.5*info.h};
        double max_pos [3] = {info.origin_x + (0.5+FluidBlock::sizeX-1)*info.h,
                              info.origin_y + (0.5+FluidBlock::sizeY-1)*info.h,
                              info.origin_z + (0.5+FluidBlock::sizeZ-1)*info.h};
        if(vSegments[s].isIntersectingWithAABB(min_pos, max_pos))
        {
          if (!hasSegments)
          {
            hasSegments = true;
            OtherSegments.resize(OtherSegments.size()+1);
          }
          OtherSegments.back().push_back(s);
        }       
      }
    }
    #pragma omp parallel
    {
      PutFishOnBlocks putfish(myFish, position, quaternion);
      #pragma omp for
      for (size_t j=0 ; j < recv_blocks[r].size(); j++)
      {
        std::vector<VolumeSegment_OBB*> S;
        for (size_t k = 0 ; k < OtherSegments[j].size() ; k++)
        {
          VolumeSegment_OBB*const ptr  = & vSegments[OtherSegments[j][k]];
            S.push_back(ptr);
        }
        if(S.size() > 0)
        {
          ObstacleBlock block;
          block.clear();
          putfish(recv_blocks[r][j].h,
                  recv_blocks[r][j].origin_x,
                  recv_blocks[r][j].origin_y,
                  recv_blocks[r][j].origin_z, &block, S);

          int kounter = 0;
          for (int iz = 0 ; iz < sizeZ ; iz++)
          for (int iy = 0 ; iy < sizeY ; iy++)
          for (int ix = 0 ; ix < sizeX ; ix++)
          {
            recv_obstacles[r][j].d[                    kounter] = block.udef[iz][iy][ix][0];
            recv_obstacles[r][j].d[sizeZ*sizeY*sizeX  +kounter] = block.udef[iz][iy][ix][1];
            recv_obstacles[r][j].d[sizeZ*sizeY*sizeX*2+kounter] = block.udef[iz][iy][ix][2];
            recv_obstacles[r][j].i[                    kounter] = block.sectionMarker[iz][iy][ix];
            kounter ++;
          }
          kounter = 0;
          for (int iz = 0 ; iz < sizeZ+2 ; iz++)
          for (int iy = 0 ; iy < sizeY+2 ; iy++)
          for (int ix = 0 ; ix < sizeX+2 ; ix++)
          {
            recv_obstacles[r][j].d[sizeZ*sizeY*sizeX*3 + kounter] = block.sdfLab[iz][iy][ix];
            kounter ++;
          }
        }
      }
    }
  }

  //Send and receive data (yes, we receive send_obstacles and send recv_obstacles)
  std::vector<MPI_Request> recv_request_obs;
  for (int r = 0 ; r < size ; r ++) if (send_obstacles[r].size() != 0)
  {
     MPI_Request req;
     recv_request_obs.push_back(req);
     MPI_Irecv(send_obstacles[r].data(), send_obstacles[r].size(), MPI_OBSTACLE, r, r*size+rank, comm, &recv_request_obs.back());
  }
  std::vector<MPI_Request> send_request_obs;
  for (int r = 0 ; r < size ; r ++) if (recv_obstacles[r].size() != 0)
  {
     MPI_Request req;
     send_request_obs.push_back(req);
     MPI_Isend(recv_obstacles[r].data(), recv_obstacles[r].size(), MPI_OBSTACLE, r, r +rank*size, comm, &send_request_obs.back());
  }

  //Compute my own blocks (that I did not send), while waiting for communication
  #pragma omp parallel
  {
    PutFishOnBlocks putfish(myFish, position, quaternion);
    #pragma omp for
    for(size_t j=counter; j< MyBlockIDs.size(); ++j)
    {
      std::vector<VolumeSegment_OBB*> S;
      for (size_t k = 0 ; k < MySegments[j].size() ; k++)
      {
        VolumeSegment_OBB*const ptr  = & vSegments[MySegments[j][k]];
          S.push_back(ptr);
      }
      if(S.size() > 0)
      {
        ObstacleBlock*const block = obstacleBlocks[MyBlockIDs[j].blockID];
        putfish(MyBlockIDs[j].h,
                MyBlockIDs[j].origin_x,
                MyBlockIDs[j].origin_y,
                MyBlockIDs[j].origin_z, block, S);
      }
    }
  }

  MPI_Waitall(send_request_obs.size(), send_request_obs.data() , MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_request_obs.size(), recv_request_obs.data() , MPI_STATUSES_IGNORE);

  counter = 0;
  for (int r = 0 ; r < size ; r ++) if (send_obstacles[r].size() != 0)
  {
    for (size_t i = 0 ; i < send_blocks[r].size() ; i ++)
    {
      ObstacleBlock*const block = obstacleBlocks[MyBlockIDs[counter + i].blockID];
      int kounter = 0;
      for (int iz = 0 ; iz < sizeZ ; iz++)
      for (int iy = 0 ; iy < sizeY ; iy++)
      for (int ix = 0 ; ix < sizeX ; ix++)
      {
        block->udef[iz][iy][ix][0]       = send_obstacles[r][i].d[                    kounter];
        block->udef[iz][iy][ix][1]       = send_obstacles[r][i].d[sizeZ*sizeY*sizeX  +kounter];
        block->udef[iz][iy][ix][2]       = send_obstacles[r][i].d[sizeZ*sizeY*sizeX*2+kounter];
        block->sectionMarker[iz][iy][ix] = send_obstacles[r][i].i[                    kounter];
        kounter ++;
      }
      kounter = 0;
      for (int iz = 0 ; iz < sizeZ+2 ; iz++)
      for (int iy = 0 ; iy < sizeY+2 ; iy++)
      for (int ix = 0 ; ix < sizeX+2 ; ix++)
      {
        block->sdfLab[iz][iy][ix] = send_obstacles[r][i].d[sizeZ*sizeY*sizeX*3 + kounter];
        kounter ++;
      }
    }
    counter += send_blocks[r].size();
  }

  #endif //load-balancing
}

void Fish::create()
{
  // STRATEGY
  // we need some things already
  // - the internal angle at the previous timestep, obtained from integrating the actual def velocities
  // (not the imposed deformation velocies, because they dont have zero ang mom)
  // - the internal angular velocity at previous timestep

  // 1. create midline
  // 2. integrate to find CoM, angular velocity, etc
  // 3. shift midline to CoM frame: zero internal linear momentum and angular momentum

  // 4. split the fish into segments (according to s)
  // 5. rotate the segments to computational frame (comp CoM and angle)
  // 6. for each Block in the domain, find those segments that intersect it
  // 7. for each of those blocks, allocate an ObstacleBlock

  // 8. put the 3D shape on the grid: SDF-P2M for sdf, normal P2M for udef
  //apply_pid_corrections();

  // 1.
  myFish->computeMidline(sim.time, sim.dt);

  // 2. & 3.
  integrateMidline();

  //CAREFUL: this func assumes everything is already centered around CM to start with, which is true (see steps 2. & 3. ...) for rX, rY: they are zero at CM, negative before and + after

  // 4. & 5.
  std::vector<VolumeSegment_OBB> vSegments = prepare_vSegments();

  // 6. & 7.
  const intersect_t segmPerBlock = prepare_segPerBlock(vSegments);

  // 8.
  writeSDFOnBlocks(vSegments);
}

void Fish::finalize() {}

void Fish::update()
{
  Obstacle::update();
}

void Fish::save(std::string filename)
{
    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);
    savestream.open(filename + ".txt");

    savestream<<sim.time<<"\t"<<sim.dt<<std::endl;
    savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
    savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
    savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
    savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
    savestream.close();
}

void Fish::restart(std::string filename)
{
  std::ifstream restartstream;
  restartstream.open(filename+".txt");
  if(!restartstream.good()){
    printf("Could not restart from file\n");
    return;
  }
  Real restart_time, restart_dt;
  restartstream >> restart_time >> restart_dt;
  restartstream >> position[0] >> position[1] >> position[2];
  restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
  restartstream >> transVel[0] >> transVel[1] >> transVel[2];
  restartstream >> angVel[0] >> angVel[1] >> angVel[2];
  restartstream.close();

  std::cout<<"RESTARTED FISH: "<<std::endl;
  std::cout<<"TIME, DT: "<<restart_time<<" "<<restart_dt<<std::endl;
  std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
  std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]
           <<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
  std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
  std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
}

CubismUP_3D_NAMESPACE_END
