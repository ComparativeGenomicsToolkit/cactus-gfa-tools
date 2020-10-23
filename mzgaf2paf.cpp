/**
 * mzgaf2paf.hpp: Make base level pairwise alignemnts from minigraph --write-mz output with the object of using them as
 *                anchors for other graph methods
 */

#include "mzgaf2paf.hpp"
#include "mzgaf.hpp"

using namespace gafkluge;
using namespace std;

void mzgaf2paf(const MzGafRecord& gaf_record, ostream& paf_stream) {

    paf_stream << gaf_record.query_name << "\t"
               << gaf_record.query_length << "\t"
               << gaf_record.query_start << "\t"
               << gaf_record.query_end << "\t"
               << (gaf_record.is_reverse ? "+" : "-") << "\t"
               << gaf_record.target_name << "\t"
               << gaf_record.target_length << "\t"
               << gaf_record.target_start << "\t"
               << gaf_record.target_end << "\t"
               << "*" << "\t"
               << "*" << "\t"
               << "*" << "\t";

    paf_stream << endl;
}

