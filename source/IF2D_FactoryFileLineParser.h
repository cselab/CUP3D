//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Wim van Rees on 04/10/14.
//

#ifndef __IF2D_ROCKS__IF2D_FactoryFileLineParser__
#define __IF2D_ROCKS__IF2D_FactoryFileLineParser__

#include "Cubism/ArgumentParser.h"

#include <algorithm>
#include <locale>
#include <sstream>
#include <utility>

class IF2D_FactoryFileLineParser: public ArgumentParser
{
protected:
    // from stackoverflow

    // trim from start
    inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    // trim from end
    inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    // trim from both ends
    inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
    }

public:

    IF2D_FactoryFileLineParser(std::istringstream & is_line)
    : ArgumentParser(0, NULL, '#') // last char is comment leader
    {
        std::string key,value;
        while( std::getline(is_line, key, '=') )
        {
            if( std::getline(is_line, value, ' ') )
            {
                // add "-" because then we can use the same code for parsing factory as command lines
                //mapArguments["-"+trim(key)] = Value(trim(value));
                mapArguments[trim(key)] = Value(trim(value));
            }
        }

        mute();
    }
};

#endif /* defined(__IF2D_ROCKS__IF2D_FactoryFileLineParser__) */
