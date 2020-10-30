/**
 * mzgaf2paf.hpp: Make base level pairwise alignemnts from minigraph --write-mz output with the object of using them as
 *                anchors for other graph methods
 */


#pragma once
#include "mzgaf.hpp"
#include <iostream>

/**
 * Convert mzgaf to paf
 */
void mzgaf2paf(const gafkluge::MzGafRecord& gaf_record,
               const gafkluge::GafRecord& parent_record,
               ostream& paf_stream,
               size_t min_gap,
               const string& target_prefix = "");
