//
//  IF2D_ObstacleFactory.cpp
//  IF2D_ROCKS
//
//  Created by Wim van Rees on 06/10/14.
//
//

#include <iostream>
#include <fstream>
#include "IF2D_ObstacleFactory.h"
#include "IF2D_FactoryFileLineParser.h"

//#include "IF2D_CylinderObstacleOperator.h"
#include "IF3D_MovingCylinderObstacleOperator.h"
#include "IF3D_ForcedCylinderObstacleOperator.h"
//#include "IF3D_StefanFishOperator.h"
#include "IF3D_CarlingFishOperator.h"
/*
#include "IF2D_EllipseObstacleOperator.h"
#include "IF2D_CylinderPairOperator.h"
 */

std::vector<IF2D_ObstacleOperator*> IF2D_ObstacleFactory::create(ArgumentParser & parser)
{
    std::vector<IF2D_ObstacleOperator*> retval;

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
    
    std::cout << "-------------   OBSTACLE FACTORY : START (" << factoryLines.size() <<" objects)   ------------" << std::endl;
    
    int k(0);
    for(auto & object : factoryLines)
    {
        std::string objectName = object.first;

        	if( objectName == "IF2D_MovingCylinder" )
        {
            retval.push_back(new IF2D_MovingCylinderObstacleOperator(grid,object.second));
        }
        else if( objectName == "IF2D_ForcedCylinder" )
        {
            retval.push_back(new IF2D_ForcedCylinderObstacleOperator(grid,object.second));
        }
        /*
        else if( objectName == "IF2D_StefanFish" )
        {
            retval.push_back(new IF2D_StefanFishOperator(grid,object.second));
        }
        */
        else if( objectName == "IF2D_CarlingFish" )
        {
            retval.push_back(new IF2D_CarlingFishOperator(grid,object.second));
        }
        /*
        else if( objectName == "IF2D_Ellipse" )
        {
            retval.push_back(new IF2D_EllipseObstacleOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF2D_CylinderPair" )
        {
            retval.push_back(new IF2D_CylinderPairOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF2D_DeadFish" )
        {
            retval.push_back(new IF2D_DeadFishOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF2D_StefanTurn" )
        {
            retval.push_back(new IF2D_StefanTurnOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF2D_StefanLearnTurn" )
        {
            retval.push_back(new IF2D_StefanLearnTurnOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF2D_Cstart" )
        {
            retval.push_back(new IF2D_CstartOperator(grid,object.second,max_bpd));
        }
        else if( objectName == "IF2D_VortexGenerator" )
        {
            retval.push_back(new IF2D_GenerateVortexOperator(grid,object.second,max_bpd));
        }
         */
        else 
        {
            std::cout << "Case " << objectName << " is not defined : skipping" << std::endl;
        }
        
        retval.back()->obstacleID = k++;
    }
    
    std::cout << "-------------   OBSTACLE FACTORY : END   ------------" << std::endl;
    
    return retval;
}
