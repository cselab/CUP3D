/*
 *  Profiler.cpp
 *  Cubism
 *
 *  Created by Diego Rossinelli on 9/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include <sys/time.h>

#include "Profiler.h"

//#include <tbb/tick_count.h>
//using namespace tbb;

void ProfileAgent::_getTime(ClockTime& time)


double ProfileAgent::_getElapsedTime(const ClockTime& tS, const ClockTime& tE)
