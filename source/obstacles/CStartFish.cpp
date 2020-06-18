//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ioannis Mandralis (ioannima@ethz.ch) adapted from code by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#include "CStartFish.h"
#include "FishLibrary.h"
#include "FishShapes.h"

#include <Cubism/ArgumentParser.h>

#include <array>
#include <cmath>
#include <sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

class ControlledCurvatureDefinedFishMidlineData : public FishMidlineData
{
public:
    // Last baseline curvature
    Real lastB3 = 0;
    Real lastB4 = 0;
    Real lastB5 = 0;
    // Last undulatory curvature
    Real lastK3 = 0;
    Real lastK4 = 0;
    Real lastK5 = 0;
    // Last tail phase
    Real lastTau = 0;
    // Last phi
    Real lastPhiUndulatory = 0;
    // Last alpha
    Real lastAlpha = 0;
    // Older baseline curvature
    Real oldrB3 = 0;
    Real oldrB4 = 0;
    Real oldrB5 = 0;
    // Older undulatory curvature
    Real oldrK3 = 0;
    Real oldrK4 = 0;
    Real oldrK5 = 0;
    // Older tail phase
    Real oldrTau = 0;
    // Older alpha
    Real oldrAlpha = 0;
    // Older phi
    Real oldrPhiUndulatory = 0;

    // Distance travelled at Tprop
    Real dTprop = 0.0;

    // first action
    bool firstAction = true;

    // Time for next action
    double t_next = 0.0;

    // Target location
    double target[2] = {0.0, 0.0};

    // Virtual origin
    double virtualOrigin[2] = {0.5, 0.5};

    // Energy expended
    double energyExpended = 0.0;
    double energyBudget = 0.0;

    // Dump time
    Real nextDump = 0.0;

    // act bools
    bool act1 = true;
    bool act2 = true;
    bool act3 = true;
    bool act4 = true;
    bool act5 = true;

    // Schedulers
    Schedulers::ParameterSchedulerVector<6> baselineCurvatureScheduler;
    Schedulers::ParameterSchedulerVector<6> undulatoryCurvatureScheduler;
    Schedulers::ParameterSchedulerScalar tauTailScheduler;
    Schedulers::ParameterSchedulerScalar phiScheduler;

protected:
    // Current curvature and curvature velocity
    Real * const rK;
    Real * const vK;
    // Current baseline curvature and curvature velocity
    Real * const rBC;
    Real * const vBC;
    // Current undulatory curvature and curvature velocity
    Real * const rUC;
    Real * const vUC;
    // Current tail phase and phase velocity
    Real tauTail;
    Real vTauTail;
    // Current phi and phi velocity
    Real phiUndulatory;
    Real vPhiUndulatory;
    // Current alpha
    Real alpha;

public:
    ControlledCurvatureDefinedFishMidlineData(Real L, Real T, Real phi, Real _h, Real _ampFac)
            : FishMidlineData(L, T, phi, _h, _ampFac), rK(_alloc(Nm)), vK(_alloc(Nm)),
              rBC(_alloc(Nm)),vBC(_alloc(Nm)), rUC(_alloc(Nm)), vUC(_alloc(Nm)),
              tauTail(0.0), vTauTail(0.0), phiUndulatory(0.0), vPhiUndulatory(0.0), alpha(0.0) {
    }

    void schedule(const Real t_current, const std::vector<double>&a)
    {
        printf("[schedule]\n");
        // Current time must be later than time at which action should be performed.
        assert(t_current >= t_rlAction);

        // Store last action into the older action placeholder
        oldrB3 = lastB3;
        oldrB4 = lastB4;
        oldrB5 = lastB5;
        oldrK3 = lastK3;
        oldrK4 = lastK4;
        oldrK5 = lastK5;
        oldrTau = lastTau;
        oldrAlpha = lastAlpha;
        oldrPhiUndulatory = lastPhiUndulatory;

        // Store the new action into the last action placeholder
        lastB3 = a[0];
        lastB4 = a[1];
        lastB5 = a[2];
        lastK3 = a[3];
        lastK4 = a[4];
        lastK5 = a[5];
        lastTau = a[6];
        lastAlpha = a[7];
        lastPhiUndulatory = a[8];

        // RL agent should output normalized curvature values as actions.
        double curvatureFactor = 1.0 / this->length;

        // Define the agent-prescribed curvature values
        const std::array<Real ,6> baselineCurvatureValues = {
                (Real)0.0 * curvatureFactor, (Real)0.0 * curvatureFactor, (Real)lastB3 * curvatureFactor,
                (Real)lastB4 * curvatureFactor, (Real)lastB5 * curvatureFactor, (Real)0.0 * curvatureFactor
        };
        const std::array<Real ,6> undulatoryCurvatureValues = {
                (Real)0.0 * curvatureFactor, (Real)0.0 * curvatureFactor, (Real)lastK3 * curvatureFactor,
                (Real)lastK4  * curvatureFactor, (Real)lastK5 * curvatureFactor, (Real)0.0 * curvatureFactor
        };

        // Using the agent-prescribed action duration get the final time of the prescribed action
        const Real actionDuration = (1 - lastAlpha) * 0.5 * this->Tperiod + lastAlpha * this->Tperiod;
        this->t_next = t_current + actionDuration;

        // Decide whether to use the current derivative for the cubic interpolation
        const bool useCurrentDerivative = true;

        // Act by scheduling a transition at the current time.
        baselineCurvatureScheduler.transition(t_current, t_current, this->t_next, baselineCurvatureValues, useCurrentDerivative);
        undulatoryCurvatureScheduler.transition(t_current, t_current, this->t_next, undulatoryCurvatureValues, useCurrentDerivative);
        tauTailScheduler.transition(t_current, t_current, this->t_next, lastTau, useCurrentDerivative);

        if (firstAction) {
            printf("FIRST ACTION %f\n", lastPhiUndulatory);
            phiScheduler.transition(t_current, t_current, this->t_next, lastPhiUndulatory, lastPhiUndulatory);
            firstAction = false;
        } else {
            printf("Next action %f\n", lastPhiUndulatory);
            phiScheduler.transition(t_current, t_current, this->t_next, lastPhiUndulatory, useCurrentDerivative);
        }

        printf("Action duration is: %f\n", actionDuration);
        printf("t_next is: %f/n", this->t_next);
        printf("Scheduled a transition between %f and %f to baseline curvatures %f, %f, %f\n", t_current, t_next, lastB3, lastB4, lastB5);
        printf("Scheduled a transition between %f and %f to undulatory curvatures %f, %f, %f\n", t_current, t_next, lastK3, lastK4, lastK5);
        printf("Scheduled a transition between %f and %f to tau %f and phi %f\n", t_current, t_next, lastTau, lastPhiUndulatory);
    }

    ~ControlledCurvatureDefinedFishMidlineData() override {
        _dealloc(rBC); _dealloc(vBC); _dealloc(rUC); _dealloc(vUC);
        _dealloc(rK); _dealloc(vK);
    }

    void computeMidline(const Real time, const Real dt) override;
};

void ControlledCurvatureDefinedFishMidlineData::computeMidline(const double t, const double dt)
{
    // Curvature control points along midline of fish, as in Gazzola et. al.
    const std::array<Real ,6> curvaturePoints = { (Real)0, (Real).2*length,
                                                  (Real).5*length, (Real).75*length, (Real).95*length, length};

    // Reproduces the 3D C-start (Gazzola et. al.)
    if (t>=0.0 && act1){
        std::vector<double> a{-1.96, -0.46, -0.56, -6.17, -3.71, -1.09, 0.65, 0.4, 0.1321};
        this->schedule(t, a);
        act1=false;
    }
    if (t>=0.7* this->Tperiod && act2){
        std::vector<double> a{0, 0, 0, -6.17, -3.71, -1.09, 0.65, 1, 0.1321};
        this->schedule(t, a);
        act2=false;
    }

//    // Reproduces the 7.30mJ energy escape
//    if (t>=0.0 && act1){
//        std::vector<double> a{-4.623, -3.75, -2.034, -1.138, -0.948, -1.374, 0.521658, 0.1651, 0.544885};
//        this->schedule(t, a);
//        act1=false;
//    }
//    if (t>=0.58255* this->Tperiod && act2){
//        std::vector<double> a{-3.082, -3.004, -1.725, -4.696, -2.979, -1.974, 0.23622, 0.1071, 0.756351};
//        this->schedule(t, a);
//        act2=false;
//    }
//    if (t>=(0.58255 + 0.553544) * this->Tperiod && act3){
//        std::vector<double> a{-1.6024, -0.9016, -2.397, -1.356, -1.633, -4.0767, 0.6017, 0.3174, 0.390727};
//        this->schedule(t, a);
//        act3=false;
//    }
//    if (t>=(0.58255 + 0.553544 + 0.658692) * this->Tperiod && act4){
//        std::vector<double> a{-1.258, -0.928, -2.5133, -3.56, -2.574, -2.9287, 0.520897, 0.2516, 0.602385};
//        this->schedule(t, a);
//        act4=false;
//    }
//    if (t>=(0.58255 + 0.553544 + 0.658692 + 0.6358) * this->Tperiod && act5){
//        std::vector<double> a{-3.04523, -2.983, -2.784, -3.868, -2.648, -2.894, 0.493, 0.3608, 0.481728};
//        this->schedule(t, a);
//        act5=false;
//    }

    // Write values to placeholders
    baselineCurvatureScheduler.gimmeValues(t, curvaturePoints, Nm, rS, rBC, vBC); // writes to rBC, vBC
    undulatoryCurvatureScheduler.gimmeValues(t, curvaturePoints, Nm, rS, rUC, vUC); // writes to rUC, vUC
    tauTailScheduler.gimmeValues(t, tauTail, vTauTail); // writes to tauTail and vTauTail
    phiScheduler.gimmeValues(t, phiUndulatory, vPhiUndulatory);

    const double curvMax = 2*M_PI/length;
#pragma omp parallel for schedule(static)
    for(int i=0; i<Nm; ++i) {

        const Real tauS = tauTail * rS[i] / length;
        const Real vTauS = vTauTail * rS[i] / length;
        const Real arg = 2 * M_PI * (t/Tperiod - tauS) + 2 * M_PI * phiUndulatory;
        const Real vArg = 2 * M_PI / Tperiod - 2 * M_PI * vTauS + 2 * M_PI * vPhiUndulatory;

        const Real curvCmd = rBC[i] + rUC[i] * std::sin(arg);
        const Real curvCmdVel = vBC[i] + rUC[i] * vArg * std::cos(arg) + vUC[i] * std::sin(arg);

        if (std::abs(curvCmd) >= curvMax) {
            rK[i] = curvCmd>0 ? curvMax : -curvMax;
            vK[i] = 0;
        } else {
            rK[i] = curvCmd;
            vK[i] = curvCmdVel;
        }

        assert(not std::isnan(rK[i]));
        assert(not std::isinf(rK[i]));
        assert(not std::isnan(vK[i]));
        assert(not std::isinf(vK[i]));
    }

    // solve Frenet to compute midline parameters
    Frenet2D::solve(Nm, rS, rK,vK, rX,rY, vX,vY, norX,norY, vNorX,vNorY);
}

CStartFish::CStartFish(SimulationData & s, ArgumentParser&p) : Fish(s, p)
{
  const double ampFac = p("-amplitudeFactor").asDouble(1.0);
  myFish = new ControlledCurvatureDefinedFishMidlineData(length, Tperiod, phaseShift, sim.maxH(), ampFac);

  std::string heightName = p("-heightProfile").asString("baseline");
  std::string  widthName = p("-widthProfile").asString("baseline");
  MidlineShapes::computeWidthsHeights(heightName, widthName, length,
                                      myFish->rS, myFish->height, myFish->width,
                                      myFish->Nm, sim.rank);

  origC[0] = position[0];
  origC[1] = position[1];
  origAng = _2Dangle;

  if(!sim.rank) printf("ControlledCurvatureDefinedFish %d %f %f %f\n",myFish->Nm, length, Tperiod, phaseShift);

}

void CStartFish::create()
{
  auto * const cFish = dynamic_cast<ControlledCurvatureDefinedFishMidlineData*>( myFish );
  if(cFish == nullptr) { printf("Someone touched my fish\n"); abort(); }
  Fish::create();
}

void CStartFish::save(std::string filename)
{
    auto * const cFish = dynamic_cast<ControlledCurvatureDefinedFishMidlineData*>( myFish );
    if(cFish == nullptr) { printf("Someone touched my fish\n"); abort(); }

    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);
    savestream.open(filename + ".txt");

    savestream<<sim.time<<"\t"<<sim.dt<<std::endl;
    savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
    savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
    savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
    savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
    savestream<<theta_internal<<"\t"<<angvel_internal<<std::endl;
    savestream.close();

    cFish->baselineCurvatureScheduler.save(filename+ "_baselineCurv");
    cFish->undulatoryCurvatureScheduler.save(filename+ "_undulatoryCurv");
    cFish->tauTailScheduler.save(filename+ "_tauTail");
    cFish->phiScheduler.save(filename+ "_phi");
}

void CStartFish::restart(std::string filename)
{
    auto * const cFish = dynamic_cast<ControlledCurvatureDefinedFishMidlineData*>( myFish );
    if(cFish == nullptr) { printf("Someone touched my fish\n"); abort(); }

    double restarted_time, restarted_dt; // timeshift, time0, l_Tp,
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    if(!restartstream.good()){
        printf("Could not restart from file\n");
        return;
    }
    restartstream >> restarted_time >> restarted_dt;
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream >> theta_internal >> angvel_internal; //  >> adjTh
    restartstream.close();

    cFish->baselineCurvatureScheduler.restart(filename+ "_baselineCurv");
    cFish->undulatoryCurvatureScheduler.restart(filename+ "_undulatoryCurv");
    cFish->tauTailScheduler.restart(filename+ "_tauTail");
    cFish->phiScheduler.restart(filename+ "_phi");

    if(!sim.rank)
    {
        std::cout<<"RESTARTED FISH: "<<std::endl;
        std::cout<<"TIME, DT: "<<restarted_time<<" "<<restarted_dt<<std::endl;
        std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
        std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
        std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
        std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
        std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
        std::cout<<"2D angle: "<<_2Dangle<<std::endl;
    }
}

CubismUP_3D_NAMESPACE_END
