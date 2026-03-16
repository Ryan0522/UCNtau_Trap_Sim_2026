#pragma once
#include "ucntrap/experiment/tracker.hpp"
#include <fstream>
#include <string>

namespace ucntrap {
 
class ResultWriter {
public:
    virtual ~ResultWriter() = default;
    virtual void write(const Result& result) = 0;
    virtual void flush() = 0;
};

class CsvResultWriter final : public ResultWriter {
public:
    explicit CsvResultWriter(const std::string& path);

    void write(const Result& result) override;
    void flush() override;

private:
    void write_header();
    
    std::ofstream out_;
    bool header_written_ = false;
};

} // namespace ucntrap 