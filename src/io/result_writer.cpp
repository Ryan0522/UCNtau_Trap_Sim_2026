#include "ucntrap/io/result_writer.hpp"

#include <fstream>
#include <memory>
#include <stdexcept>
#include <iomanip>

namespace ucntrap {

CsvResultWriter::CsvResultWriter(const std::string& path, int precision)
    : out_(path), buffer_(1 << 20) {
    if (!out_.is_open()) {
        throw std::runtime_error(
            "CsvResultWriter: unable to open out file '" + path + "'"
        );
    }

    out_ << std::scientific << std::setprecision(precision);
    out_.rdbuf()->pubsetbuf(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
}

void CsvResultWriter::write_header() {
    if (header_written_) {
        return;
    }

    out_ << "t_start,t_final,e_start,e_final,"
         << "x_start,y_start,z_start,vx_start,vy_start,vz_start,"
         << "x_final,y_final,z_final,vx_final,vy_final,vz_final,"
         << "z_off,n_hit,n_hit_house_low,n_hit_house_high,death_time,code\n";

    if (!out_) {
        throw std::runtime_error("CsvResultWriter: failed while writing header");
    } 

    header_written_ = true;
}

void CsvResultWriter::write(const Result& result) {
    write_header();

    out_ << result.t_start << ','
         << result.t_final << ','
         << result.e_start << ','
         << result.e_final << ','
         << result.x_start << ','
         << result.y_start << ','
         << result.z_start << ','
         << result.vx_start << ','
         << result.vy_start << ','
         << result.vz_start << ','
         << result.x_final << ','
         << result.y_final << ','
         << result.z_final << ','
         << result.vx_final << ','
         << result.vy_final << ','
         << result.vz_final << ','
         << result.z_off << ','
         << result.n_hit << ','
         << result.n_hit_house_low << ','
         << result.n_hit_house_high << ','
         << result.death_time << ','
         << result.code
         << '\n';

    if (!out_) {
        throw std::runtime_error("CsvResultWriter: failed while writing row");
    }
}

void CsvResultWriter::flush() {
    out_.flush();
    if (!out_) {
        throw std::runtime_error("CsvResultWriter: flush failed");
    }
}

} // namespace ucntrap 