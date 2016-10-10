//
//  IF3D_ObstacleFactory.cpp
//  IF3D_ROCKS
//
//  Created by Wim van Rees on 06/10/14.
//
//

#include <iostream>
#include <fstream>
#include "IF3D_ObstacleFactory.h"
#include "IF2D_FactoryFileLineParser.h"

//#include "IF3D_CylinderObstacleOperator.h"
#include "IF3D_SphereObstacleOperator.h"
#include "IF3D_ForcedSphereObstacleOperator.h"
#include "IF3D_StefanFishOperator.h"
#include "IF3D_DeadFishOperator.h"
#include "IF3D_CarlingFishOperator.h"
/*
#include "IF3D_EllipseObstacleOperator.h"
#include "IF3D_CylinderPairOperator.h"
 */

std::vector<IF3D_ObstacleOperator*> IF3D_ObstacleFactory::create(ArgumentParser & parser)
{
    std::vector<IF3D_ObstacleOperator*> retval;

    // Read parser information
    parser.unset_strict_mode();
    const std::string factoryFileName = parser("-factory").asString("factory");
    assert(factoryFileName != "");

    // open the file
    std::ifstream is_file(factoryFileName.c_str());

    // here we store the data per object
    std::vector<std::pair<std::string,IF2D_FactoryFileLineParser>> factoryLines;
    std::string line;

    while( std::getline(is_file, line) ) {
        std::istringstream is_line(line);
        std::string ID;
        is_line >> ID;
        if(ID.empty()) continue;
        IF2D_FactoryFileLineParser ffparser(is_line);
        factoryLines.push_back(make_pair(ID,ffparser));
    }
    if(rank==0)
    std::cout << "-------------   OBSTACLE FACTORY : START (" << factoryLines.size() <<" objects)   ------------" << std::endl;
    
    int k(0);
    for(auto & object : factoryLines)
    {
        std::string objectName = object.first;

        	if( objectName == "IF3D_Sphere" )
        {
            retval.push_back(new IF3D_SphereObstacleOperator(grid,object.second));
        }
        else if( objectName == "IF3D_ForcedSphere" )
        {
            retval.push_back(new IF3D_ForcedSphereObstacleOperator(grid,object.second));
        }
        else if( objectName == "IF3D_DeadFish" )
        {
            retval.push_back(new IF3D_DeadFishOperator(grid,object.second));
        }
        else if( objectName == "IF3D_StefanFish" )
        {
            retval.push_back(new IF3D_StefanFishOperator(grid,object.second));
        }
        else if( objectName == "IF3D_CarlingFish" )
        {
            retval.push_back(new IF3D_CarlingFishOperator(grid,object.second));
        }
        /*
        else if( objectName == "IF3D_Ellipse" )
        {
            retval.push_back(new IF3D_EllipseObstacleOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF3D_CylinderPair" )
        {
            retval.push_back(new IF3D_CylinderPairOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF3D_DeadFish" )
        {
            retval.push_back(new IF3D_DeadFishOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF3D_StefanTurn" )
        {
            retval.push_back(new IF3D_StefanTurnOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF3D_StefanLearnTurn" )
        {
            retval.push_back(new IF3D_StefanLearnTurnOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF3D_Cstart" )
        {
            retval.push_back(new IF3D_CstartOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF3D_VortexGenerator" )
        {
            retval.push_back(new IF3D_GenerateVortexOperator(grid,object.second,max_bpd));
        }
         */
        else 
        {
            if(rank==0)
            std::cout << "Case " << objectName << " is not defined : skipping" << std::endl;
        }
        
        retval.back()->obstacleID = k++;
    }

    if(rank==0)
    std::cout << "-------------   OBSTACLE FACTORY : END   ------------" << std::endl;
    
    return retval;
}
