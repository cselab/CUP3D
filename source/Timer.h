/*
 *  Timer.h
 *
 *	Timer class using C++11 chrono
 *
 *  Created by Christian Conti on 11/12/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <chrono>

class Timer
{
    std::chrono::time_point<std::chrono::high_resolution_clock> t_start, t_end;
	
public:
	
	void start()
	{
		t_start = std::chrono::high_resolution_clock::now();
	}
	
	double stop()
	{
		t_end = std::chrono::high_resolution_clock::now();
		
		return std::chrono::duration<double>(t_end-t_start).count();
	}
};
