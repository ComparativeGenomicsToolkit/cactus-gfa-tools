// Taken from https://raw.githubusercontent.com/edawson/tinyGFA/master/tinygfa.hpp
// https://github.com/edawson/tinyGFA/commit/e6828ec888591e65b4f4c9a58c601043111bb057

#ifndef TGFA_H_DEF
#define TGFA_H_DEF
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdint>
#include <functional>

#include <pliib.hpp>


namespace tgfa{



    enum GFA_LINE_TYPES {HEADER_LINE,SEQUENCE_LINE,FRAGMENT_LINE,EDGE_LINE,GAP_LINE,GROUP_LINE,PATH_LINE,LINK_LINE,CONTAINED_LINE,WALK_LINE};

    static inline int determine_line_type(const char* line){
        if (line[0] == 'H'){
            return HEADER_LINE;
        }
        else if (line[0] == 'S'){
            return SEQUENCE_LINE;
        }
        else if (line[0] == 'E'){
            return EDGE_LINE;
        }
        else if (line[0] == 'L'){
            return LINK_LINE;
        }
        else if (line[0] == 'C'){
            return CONTAINED_LINE;
        }
        else if (line[0] == 'U' || line[0] == 'O'){
            return GROUP_LINE;
        }
        else if (line[0] == 'P'){
            return PATH_LINE;
        }
        else if (line[0] == 'W'){
            return WALK_LINE;
        }
        else if (line[0] == 'G'){
            return GAP_LINE;
        }
        else if (line[0] == 'F'){
            return FRAGMENT_LINE;
        }
        else{
            return -1;
        }
    };



    enum val_types {STRING_TYPE, INT_TYPE, FLOAT_TYPE};
    struct opt_elem{
        char* opt_id = nullptr;
        std::uint8_t type;
        char* val = nullptr;
        void clear(){
            if (opt_id != nullptr)
                delete [] opt_id;
            if (val != nullptr)
                delete [] val;
        }
        
        std::string get_id(){
            return std::string(opt_id);
        }
        std::string get_val(){
            return std::string(val);
        }
        std::uint8_t get_type(){
            return type;
        }
        
    };
    struct header_elem{
        char* key = nullptr;
        std::uint8_t type;
        char* val = nullptr;
    };
    struct sequence_elem{
        char* seq_id = nullptr;
        char* seq = nullptr;
        std::uint32_t seq_length;
        std::vector<opt_elem> tags;
        void clear(){
            if (seq_id != nullptr)
                delete [] seq_id;
            if (seq != nullptr)
                delete [] seq;
        }
        void set(char**& splits, std::size_t& split_count, std::size_t*& split_lens, int spec = 2){
            //clear();
            seq_id = splits[1];
            if (spec == 2){
                seq_length = std::stoull(splits[2]);
                seq = splits[3];
            }
            else if (spec == 1){
                seq = splits[2];
                seq_length = strlen(seq);
            }
            for (std::size_t i = 3; i < split_count; ++i){
            }


        }
        sequence_elem(){

        }
        sequence_elem(char**& splits, std::size_t split_count, std::size_t*& split_lens){
            set(splits, split_count, split_lens);
        }
        ~sequence_elem(){
           // clear();
        }

        inline std::ostream& output(std::ostream& os, int spec = 2){
            if (spec == 2){
                os << 'S' << '\t' << seq_id << '\t' << seq_length << '\t' << seq << '\t';
            }
            else if (spec == 1){
                os << 'S' << '\t' << seq_id  << '\t' << seq << '\t';

            }
            else {
                std::cerr << "Invalid spec: " << spec << ". Please provide a valid spec (1 or 2)" << std::endl;
                exit(9);
            }
            return os;
        }

        inline std::string to_string(int spec = 2){
            std::stringstream st;
            output(st, spec);
            return st.str();
        }
    };

    std::uint8_t EDGE_T_SINK_END_OFFSET = 1;
    std::uint8_t EDGE_T_SINK_START_OFFSET = 2;
    std::uint8_t EDGE_T_SOURCE_END_OFFSET = 4;
    std::uint8_t EDGE_T_SOURCE_START_OFFSET = 8;
    struct edge_elem{
        char* edge_id = nullptr;
        char* source_id = nullptr;
        char* sink_id = nullptr;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        std::pair<uint32_t, uint32_t> source_pos;
        std::pair<uint32_t, uint32_t> sink_pos;
        std::uint8_t ends;
        char* alignment = nullptr;
        std::vector<opt_elem> tags;

        edge_elem(){

        }
        ~edge_elem(){
        }
        edge_elem(char**& tokens, std::size_t& token_count, std::size_t*& token_sizes){
            set(tokens, token_count, token_sizes);
        }
        bool sink_end(){
            
            return true;
        };
        bool sink_begin(){
            return true;
        };
        bool source_end(){
            return true;
        };
        bool source_begin(){
            return true;
        };
        void set_end_begin(std::uint8_t& ends, std::uint8_t offset){
            ends |= 1<<offset;
        }
        void unset_end_begin(std::uint8_t& ends, std::uint8_t offset){
            ends &= ~(1<<offset);
        }
        void set(char**& tokens, std::size_t& token_count, std::size_t*& token_sizes){
            if (token_count > 0){
                // Check if the line is a GFA1 L/C line,
                // or a GFA2 E line
                int line_type = determine_line_type(tokens[0]);
                if (line_type == EDGE_LINE){
                    edge_id = tokens[1];
                    source_orientation_forward = tokens[2][token_sizes[2] - 1] == '+';
                    sink_orientation_forward = tokens[3][token_sizes[3] - 1] == '+';
                    pliib::strcopy(tokens[2], source_id);
                    pliib::strip(source_id, token_sizes[2], '+');
                    pliib::strip(source_id, strlen(source_id), '-');
                    pliib::strcopy(tokens[3], sink_id);
                    pliib::strip(sink_id, token_sizes[3], '+');
                    pliib::strip(sink_id, strlen(sink_id), '-');

                    tokens[4][token_sizes[4] - 1] == '$' ? set_end_begin(ends,EDGE_T_SOURCE_START_OFFSET) : unset_end_begin(ends,EDGE_T_SOURCE_START_OFFSET);
                    tokens[5][token_sizes[5] - 1] == '$' ? set_end_begin(ends,EDGE_T_SOURCE_END_OFFSET) : unset_end_begin(ends,EDGE_T_SOURCE_END_OFFSET);
                    tokens[6][token_sizes[6] - 1] == '$' ? set_end_begin(ends,EDGE_T_SINK_START_OFFSET) : unset_end_begin(ends,EDGE_T_SINK_START_OFFSET);
                    tokens[7][token_sizes[7] - 1] == '$' ? set_end_begin(ends,EDGE_T_SINK_END_OFFSET) : unset_end_begin(ends,EDGE_T_SINK_END_OFFSET);

                    alignment = tokens[8];
                    for (std::size_t i = 9; i < token_count; ++i){
                        opt_elem o;
                        tags.push_back(o);
                    }

                    

                }
                else if(line_type == CONTAINED_LINE){
                    std::cerr << "Containment lines not implemented [tinyfa::edge_elem]" << std::endl;
                    throw 20;
                }
                else if (line_type == LINK_LINE){
                    source_id = tokens[1];
                    source_orientation_forward = *tokens[2] == '+';
                    sink_id = tokens[3];
                    sink_orientation_forward = *tokens[4] == '+';
                    set_end_begin(ends, EDGE_T_SOURCE_START_OFFSET);
                    set_end_begin(ends, EDGE_T_SOURCE_END_OFFSET);
                    set_end_begin(ends, EDGE_T_SINK_START_OFFSET);
                    set_end_begin(ends, EDGE_T_SINK_END_OFFSET);
                }
                else{
                    std::cerr << "Invalid line type " << tokens[0] << std::endl;
                    throw 20;
                }
            }
        }
        void clear(){
            delete [] source_id;
            delete [] sink_id;
        }
    };
    struct group_elem{
        char* group_id = nullptr;
        bool ordered = true;
        std::uint64_t segment_count;
        char** segment_ids = nullptr;
        std::vector<opt_elem> tags;

        group_elem(char**& tokens, std::size_t& token_count, std::size_t*& token_sizes){
            set(tokens, token_count, token_sizes);
        }
        void set(char**& tokens, std::size_t& token_count, std::size_t*& token_sizes){
            group_id = tokens[1];
            std::size_t num_seg_splits;
            std::size_t* seg_split_lens;
            pliib::split(tokens[2], ',', segment_ids, num_seg_splits, seg_split_lens);
            for (std::size_t i = 0; i < num_seg_splits; ++i){
                pliib::slice(segment_ids[i], 0, seg_split_lens[i] - 1, segment_ids[i]);
            }
            ordered = tokens[0][0] == 'O';
            segment_count = num_seg_splits;
        }
        ~group_elem(){
            //delete [] group_id;
            //delete [] segment_ids;
        }
    };

    struct gfa_func_t{
        std::vector<std::function<void(sequence_elem)>> seq_funcs;
        std::vector<std::function<void(edge_elem)>> edge_funcs;
        std::vector<std::function<void(group_elem)>> group_funcs;
    };

    struct gfa_stat_t{
            std::uint64_t sequence_count = 0;
            std::uint64_t edge_count = 0;
            std::uint64_t group_count = 0;
            std::uint64_t path_count = 0;
            double N50 = 0.0;
            double N90 = 0.0;
            std::string to_string(){
                std::stringstream st;
                st << "Number of sequences (nodes): " << sequence_count << std::endl <<
                    "Number of edges: " << edge_count << std::endl <<
                    "Number of groups: " << group_count << std::endl <<
                    "Number of paths (ordered groups): " << path_count;
                return st.str();
            }
    };

    // Cheers to pfultz2.com/blog, 2014 09 02
    template<typename T>
    typename std::remove_reference<T>::type *addr(T &&t){
        return &t;
    }
    const constexpr auto null_func = true ? nullptr : addr([](auto x){});


    inline bool parse_gfa_file(const char* filename,
            auto& seq_func,
            auto& edge_func,
            auto& group_func,
            gfa_stat_t& stats,
            int spec = 2){

        return parse_gfa_file(filename, seq_func, true, edge_func, true, group_func, true, stats, spec);
        
    }

    inline bool parse_gfa_file(const char* filename,
            auto& seq_func,
            bool process_seqs,
            auto& edge_func,
            bool process_edges,
            auto& group_func,
            bool process_groups,
            gfa_stat_t& stats,
            int spec = 2){

        std::ifstream instream;
        instream.open(filename, std::ifstream::in);

        if (!instream.good()){
            std::cerr << "Error: file " << filename << " could not be properly opened." << std::endl;
            exit(9);
        }

        return parse_gfa_file(instream, seq_func, process_seqs, edge_func, process_edges, group_func, process_groups, stats, spec);


    }




    inline bool parse_gfa_file(std::ifstream& instream,
            auto& seq_func,
            auto& edge_func,
            auto& group_func,
            gfa_stat_t& stats,
            int spec){
            return parse_gfa_file(instream, seq_func, true, edge_func, true, group_func, true, stats, spec);
    }
    inline bool parse_gfa_file(std::ifstream& instream,
            auto& seq_func,
            bool process_seqs,
            auto& edge_func,
            bool process_edges,
            auto& group_func,
            bool process_groups,
            gfa_stat_t& stats,
            int spec){

            bool ret = true;

            if (!instream.good()){
                std::cerr << "Error: input stream failure." << std::endl;
                exit(9);
            }

            std::streamsize max_ln_size = 250000;
            char* line = new char[max_ln_size];
            char** splits;
            std::size_t split_count;
            std::size_t* split_lens;

            sequence_elem s;
            edge_elem e;
            //group_elem g;
            while(instream.getline(line, max_ln_size)){
                auto line_type = determine_line_type(line);
                if (line_type == SEQUENCE_LINE){
                    if (process_seqs){
                        pliib::split(line, '\t', splits, split_count, split_lens);
                        s.set(splits, split_count, split_lens, spec);
                        seq_func(s);
                        pliib::destroy_splits(splits, split_count, split_lens);
                        ++stats.sequence_count;
                    }
                }
                else if (line_type == EDGE_LINE || line_type == LINK_LINE || line_type == CONTAINED_LINE){
                    if (process_edges){
                        pliib::split(line, '\t', splits, split_count, split_lens);
                        e.set(splits, split_count, split_lens);
                        edge_func(e);
                        e.clear();
                        pliib::destroy_splits(splits, split_count, split_lens);
                    }
                    ++stats.edge_count;
                }
                else if (line_type == GROUP_LINE){
                    if (process_groups){
                        pliib::split(line, '\t', splits, split_count, split_lens);
                        group_elem g(splits, split_count, split_lens);
                        group_func(g);
                        pliib::destroy_splits(splits, split_count, split_lens);
                        ++stats.group_count;
                        if (g.ordered){
                            ++stats.path_count;
                        }
                    }
                }
            }
            delete [] line;

            return ret;
    }
}
#endif
