//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch) in May 2018.
//

#include "BufferedLogger.h"
#include <fstream>

namespace cubismup3d {

BufferedLogger logger;

static constexpr int AUTO_FLUSH_COUNT = 100;

void BufferedLogger::flush(BufferedLogger::container_type::iterator it) {
    std::ofstream savestream;
    savestream.open(it->first, std::ios::app | std::ios::out);
    savestream << it->second.stream.rdbuf();
    savestream.close();
    it->second.requests_since_last_flush = 0;
}

std::stringstream& BufferedLogger::get_stream(const std::string &filename) {
    auto it = files.find(filename);
    if (it != files.end()) {
        if (++it->second.requests_since_last_flush == AUTO_FLUSH_COUNT)
            flush(it);
        return it->second.stream;
    } else {
        // With request_since_last_flush == 0,
        // the first flush will have AUTO_FLUSH_COUNT frames.
        auto new_it = files.emplace(filename, Stream()).first;
        return new_it->second.stream;
    }
}

}
