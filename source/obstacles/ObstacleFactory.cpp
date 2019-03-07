//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#include <iostream>
#include <fstream>
#include "obstacles/ObstacleFactory.h"
#include "obstacles/extra/FactoryFileLineParser.h"

#include "obstacles/CarlingFish.h"
#include "obstacles/Cylinder.h"
#include "obstacles/Ellipsoid.h"
#include "obstacles/ExternalObstacle.h"
#include "obstacles/Naca.h"
#include "obstacles/Plate.h"
#include "obstacles/Sphere.h"
#include "obstacles/StefanFish.h"

//#include "obstacles/IF3D_ElasticFishOperator.h"
//#include "obstacles/IF3D_CylinderPairOperator.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

ObstacleFactory::VectorType
ObstacleFactory::create(cubism::ArgumentParser & parser)
{
  // Read parser information
  parser.unset_strict_mode();
  const std::string factory_filename = parser("-factory").asString("factory");
  const std::string factory_content = parser("-factory-content").asString("");

  assert(factory_filename != "" || factory_content != "");

  std::stringstream stream(factory_content);
  if (!factory_filename.empty()) {
      // https://stackoverflow.com/questions/132358/how-to-read-file-content-into-istringstream
      // Good enough solution.
      std::ifstream file(factory_filename);
      if (file.is_open()) {
          stream << '\n';
          stream << file.rdbuf();
      }
  }
  if (sim.rank == 0)
      printf("Factory (file) + factory (cmdline argument):\n%s\n\n", stream.str().c_str());
  // here we store the data per object
  std::vector<std::pair<std::string, FactoryFileLineParser>> factoryLines;
  std::string line;

  while (std::getline(stream, line)) {
      std::istringstream line_stream(line);
      std::string ID;
      line_stream >> ID;
      if(ID.empty() or ID[0]=='#') continue; // Comments and empty lines ignored
      FactoryFileLineParser ffparser(line_stream);
      factoryLines.push_back(make_pair(ID,ffparser));
  }
  if(factoryLines.size() == 0) {
    if(sim.rank==0)
      std::cout<<"OBSTACLE FACTORY did not create any obstacles."<<std::endl;
    return {};
  }
  if(sim.rank==0)
  std::cout << "-------------   OBSTACLE FACTORY : START (" << factoryLines.size() <<" objects)   ------------" << std::endl;

  VectorType retval;
  retval.reserve(factoryLines.size());
  int k = 0;
  for(auto & object : factoryLines)
  {
    std::string objectName = object.first;

      if( objectName == "Sphere" )
    {
      retval.emplace_back(std::make_shared<Sphere>(sim, object.second));
    }
    else if( objectName == "StefanFish" )
    {
      retval.emplace_back(std::make_shared<StefanFish>(sim, object.second));
    }
    else if( objectName == "CarlingFish" )
    {
      retval.emplace_back(std::make_shared<CarlingFish>(sim, object.second));
    }
    else if( objectName == "Naca" )
    {
      retval.emplace_back(std::make_shared<Naca>(sim, object.second));
    }
    else if( objectName == "Cylinder" )
    {
      retval.emplace_back(std::make_shared<Cylinder>(sim, object.second));
    }
    else if( objectName == "Plate" )
    {
      retval.emplace_back(std::make_shared<Plate>(sim, object.second));
    }
    else if( objectName == "ExternalObstacle" )
    {
      retval.emplace_back(std::make_shared<ExternalObstacle>(sim, object.second));
    }
    else if( objectName == "Ellipsoid" )
    {
      retval.emplace_back(std::make_shared<Ellipsoid>(sim, object.second));
    }
    /*
    else if( objectName == "ElasticFish" )
    {
      retval.emplace_back(std::make_shared<ElasticFish>(grid, object.second, Uinf));
    }
    else if( objectName == "CylinderPair" )
    {
        retval.emplace_back(std::make_shared<CylinderPair>(grid, object.second, max_bpd));
    }
    else if( objectName == "Cstart" )
    {
        retval.emplace_back(std::make_shared<Cstart>(grid, object.second, max_bpd));
    }
    else if( objectName == "VortexGenerator" )
    {
        retval.emplace_back(std::make_shared<GenerateVortex>(grid, object.second, max_bpd));
    }
     */
    else
    {
      if (sim.rank == 0) {
        std::cout<<"Case "<<objectName<<" is not defined: aborting"<<std::endl;
        abort();
      }
    }

    retval.back()->obstacleID = k++;
  }

  if(sim.rank==0)
  std::cout<<"-------------   OBSTACLE FACTORY : END   ------------"<<std::endl;

  return retval;
}

CubismUP_3D_NAMESPACE_END
