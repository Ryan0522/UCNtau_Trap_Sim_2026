#pragma once
#include "ucntrap/trace.hpp"
#include <string>

namespace ucntrap {
    
Trace load_trace(const std::string& x_file,
                 const std::string& y_file,
                 const std::string& z_file);
                 
} // namespace ucntrap 