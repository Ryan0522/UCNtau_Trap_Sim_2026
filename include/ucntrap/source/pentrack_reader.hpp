#pragma once
#include "ucntrap/config.hpp"
#include "ucntrap/source/initial_condition_source.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace ucntrap {
    
class PenTrackReader final : public InitialConditionSource {
public:
    explicit PenTrackReader(const SimulationConfig& config);
    explicit PenTrackReader(std::string filename,
                            std::size_t ntraj,
                            std::size_t skip_lines = 0);

    bool has_next() const override;
    State next() override;

    static State parse_line(const std::string& line);

private:
    std::vector<State> states_;
    std::size_t index_ = 0;
};

} // namespace ucntra 
