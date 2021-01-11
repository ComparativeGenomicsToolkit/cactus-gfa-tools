#pragma once

#include <string>

// convert paf to lastz cigar
// if use_mapq is true, then the MAPQ will be used for lastz CIGAR score, otherwise
// it will try to use the value in the AS tag (like paftools) and write 0 if it's not found
// returns lastz line paired with secondary flag
std::pair<std::string, bool> paf2lastz(const std::string& paf_line, bool use_mapq);
