//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "../ObstacleBlock.h"
#include "../SimulationData.h"

#include <array>

/*
 * HOW OBSTACLES WORK
 *
 * Each obstacle type T defines 1) an operator `T`, 2) an arguments struct
 * `TArguments` and 3) a struct `ObstacleAndTArguments` combined struct.
 *
 * 1) The operator `T` derives from a base operator `Obstacle` and from the
 *    struct `TArguments`.
 * 2) The struct `TArguments` is a stand-alone struct that defines only the
 *    additional properties of the obstacle type T.
 * 3) The struct `ObstacleAndTArguments`, used only for construction (from
 *    Python bindings), derives from `ObstacleArguments` and `TArguments`, and
 *    defined a constructor these two components.
 *
 * The operator `T` then defines the following constructor:
 *      T(SimulationData &s, const ObstacleAndTArguments &args)
 *        : Obstacle(s, args), TArguments(args) { ... }
 *
 * See `Sphere.h` for an example.
 *
 *
 * POTENTIAL ALTERNATIVE IMPLEMENTATION
 *
 * Make `TArguments` a derived class of `ObstacleArguments`.
 * With that, `Obstacle` must be a template class, e.g.
 *      `class Sphere : ObstacleImpl<Sphere, SphereArguments> { ... }`.
 *
 * Where `ObstacleImpl<>` derives from an abstract `Obstacle`.
 */


// forward declaration of derived class for visitor

namespace cubism { class ArgumentParser; }

CubismUP_3D_NAMESPACE_BEGIN

class Obstacle;
class ObstacleVector;

/*
 * Structure containing all externally configurable parameters of a base obstacle.
 */
struct ObstacleArguments
{
  Real length = 0.0;
  Real planarAngle = 0.0;
  std::array<Real, 3> position = {{0.0, 0.0, 0.0}};
  std::array<Real, 4> quaternion = {{1.0, 0.0, 0.0, 0.0}};
  std::array<Real, 3> enforcedVelocity = {{0.0, 0.0, 0.0}};  // Only if bForcedInSimFrame.
  std::array<bool, 3> bForcedInSimFrame = {{false, false, false}};
  std::array<bool, 3> bFixFrameOfRef = {{false, false, false}};
  bool bFixToPlanar = false;
  bool bBreakSymmetry = false;

  ObstacleArguments() = default;

  /* Convert human-readable format into internal representation of parameters. */
  ObstacleArguments(const SimulationData & sim, cubism::ArgumentParser &parser);
};


struct ObstacleVisitor
{
  virtual ~ObstacleVisitor() {}
  virtual void visit(Obstacle* const obstacle) = 0;
};

class Obstacle
{
protected:
  SimulationData & sim;
  FluidGridMPI * const grid = sim.grid;
  bool printedHeaderVels = false;
public:
  std::vector<ObstacleBlock*> obstacleBlocks;

  int obstacleID=0; //each obstacle has a unique ID

  //There are two frames of reference, denoted by (I) and (II):
  // (I):  An absolute frame that is not moving.
  // (II): A frame that is moving with sim.uinf. Generally, this frame follows the obstacles.
  //       sim.uinf is the average velocity of obstacles with bFixFrameOfRef = true
  //       (with the addition of a far field velocity field, if the user wants that).

  Real absPos  [3] = {0,0,0};     //position of obstacle in frame (I)
  Real position[3] = {0,0,0};     //position of obstacle in frame (II)
  Real quaternion[4] = {1,0,0,0}; //orientation (same for both frames)
  Real transVel[3] = {0,0,0};     //translational velocity in frame (I)
  Real angVel  [3] = {0,0,0};     //angular velocity
  Real mass;                      //obstacle mass
  Real length;                    //characteristic length of obstacle
  Real J[6] = {0,0,0,0,0,0};      //moments of inertia matrix (Jxx,Jyy,Jzz,Jxy,Jxz,Jyz}

  std::array<bool, 3> bFixFrameOfRef    = {{false, false, false}};//set to true if 'camera' will follow the obstacle in that direction
  std::array<bool, 3> bForcedInSimFrame = {{false, false, false}};//set to true if obstacle is forced
  std::array<bool, 3> bBlockRotation    = {{false, false, false}};//set to true if obstacle is not allowed to rotate (forced)
  Real transVel_imposed[3]= {0,0,0}; //prescribed velocity (if the obstacle is forced)

  //auxiliary arrays used for 2nd order time integration of obstacle's position
  Real old_position  [3] =   {0,0,0};
  Real old_absPos    [3] =   {0,0,0};
  Real old_quaternion[4] = {1,0,0,0};

  //The obstacle's mass, linear momentum, angular momentum, center of mass and
  //moments of inertia are computed from the chi field and stored here.
  //Then, a 6x6 linear system is solved to determine transVel and angVel.
  Real penalM;
  std::array<Real,3> penalLmom = {0,0,0};
  std::array<Real,3> penalAmom = {0,0,0};
  std::array<Real,3> penalCM   = {0,0,0};
  std::array<Real,6> penalJ    = {0,0,0,0,0,0};

  //Translational and angular velocities computed by solving the 6x6 linear system.
  //If the obstacle is not forced, transVel and angVel are set equal to transVel_computed
  //and angVel_computed. Otherwise, the forces and torques acting on the forced obstacle
  //can be computed as force = mass * (transVel_computed - transVel) / sim.dt
  std::array<Real,3> transVel_computed = {0,0,0};
  std::array<Real,3> angVel_computed   = {0,0,0};

  Real centerOfMass[3] = {0,0,0}; //center of mass (used to resolve collisions), computed from chi

  bool bBreakSymmetry = false; //if true obstacle is forced for a short period, to break symmetric flow conditions

  //Quantities of Interest (forces, drag etc.)
  std::array<Real,3> force  = {0,0,0};
  std::array<Real,3> torque = {0,0,0};
  Real gamma[3]={0,0,0};
  Real surfForce[3]={0,0,0};
  Real presForce[3]={0,0,0};
  Real viscForce[3]={0,0,0};
  Real surfTorque[3]={0,0,0};
  Real drag=0, thrust=0, Pout=0, PoutBnd=0, pLocom=0;
  Real defPower=0, defPowerBnd=0, Pthrust=0, Pdrag=0, EffPDef=0, EffPDefBnd=0;

  //Used to 'stabilize' some computation in KernelIntegrateUdefMomenta from ObstaclesCreate.cpp
  std::array<Real,3> transVel_correction={0,0,0}, angVel_correction={0,0,0};

protected:
  //functions to save quantities of interest to files
  virtual void _writeComputedVelToFile();
  virtual void _writeDiagForcesToFile();
  virtual void _writeSurfForcesToFile();

public:
  Obstacle(SimulationData& s, const ObstacleArguments &args);
  Obstacle(SimulationData& s, cubism::ArgumentParser &parser);
  Obstacle(SimulationData& s) : sim(s) {  }

  virtual void Accept(ObstacleVisitor * visitor);

  virtual void computeVelocities();//solve the 6x6 linear system to get transVel and angvel
  virtual void computeForces();    //compute quantities of interest for this obstacle
  virtual void update();           //time integration of position and orientation
  virtual void save(std::string filename = std::string());   //save information for restart
  virtual void restart(std::string filename = std::string());//read information for restart
  virtual void create();  //additional stuff to be done when creating an obstacle (optional)
  virtual void finalize();//additional stuff to be done when deleting an obstacle (optional)
  std::array<Real,3> getTranslationVelocity() const;
  std::array<Real,3> getAngularVelocity() const;
  std::array<Real,3> getCenterOfMass() const;

  std::vector<ObstacleBlock*>  getObstacleBlocks() const {return  obstacleBlocks;}
  std::vector<ObstacleBlock*>* getObstacleBlocksPtr()    {return &obstacleBlocks;}

  virtual ~Obstacle()
  {
    for(auto & entry : obstacleBlocks) {
      if(entry != nullptr) {
          delete entry;
          entry = nullptr;
      }
    }
    obstacleBlocks.clear();
  }

  // driver to execute finite difference kernels either on all points relevant
  // to the mass of the obstacle (where we have char func) or only on surface
  template<typename T>
  void create_base(const T& kernel)
  {
    for(auto & entry : obstacleBlocks) {
      if(entry == nullptr) continue;
      delete entry;
      entry = nullptr;
    }
    std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
    obstacleBlocks.resize(vInfo.size(), nullptr);

    #pragma omp parallel for schedule(dynamic, 1)
    for(size_t i=0; i<vInfo.size(); i++) {
      const cubism::BlockInfo& info = vInfo[i];
      const FluidBlock &b = *(FluidBlock *)info.ptrBlock;
      if(kernel.isTouching(info, b)) {
        assert(obstacleBlocks[info.blockID] == nullptr);
        obstacleBlocks[info.blockID] = new ObstacleBlock();
        obstacleBlocks[info.blockID]->clear(); //memset 0
        kernel(info, obstacleBlocks[info.blockID]);
      }
    }
  }
};

CubismUP_3D_NAMESPACE_END
