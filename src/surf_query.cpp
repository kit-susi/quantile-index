#include "surf/config.hpp"
#include "surf/indexes.hpp"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace sdsl;
using namespace surf;

typedef struct cmdargs {
    std::string collection_dir = "";
    std::string query_file = "";
    uint64_t k = 10;
    bool verbose = false;
    bool multi_occ = false;
    bool match_only = false;
    bool intersection = false;
    uint64_t snippet_size = 0;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout, "%s -c <collection directory> -q <query file> other options\n", program);
    fprintf(stdout, "where\n");
    fprintf(stdout, "  -c <collection directory> : the directory the collection is stored.\n");
    fprintf(stdout, "  -q <query file>   : the queries to be performed.\n");
    fprintf(stdout, "  -k <top-k>        : the top-k documents to be retrieved for each query.\n");
    fprintf(stdout, "  -v <verbose>      : verbose mode.\n");
    fprintf(stdout, "  -m <multi_occ>    : only retrieve documents which contain the term more than once.\n");
    fprintf(stdout, "  -o <only_match>   : only match pattern; no document retrieval.\n");
    fprintf(stdout, "  -s <snippet_size> : extract snippets of size snippet_size.\n");
};


int stick_this_thread_to_core(int core_id) {
    int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (core_id < 0 || core_id >= num_cores)
        return EINVAL;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    pthread_t current_thread = pthread_self();
    return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

cmdargs_t
parse_args(int argc, char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.query_file = "";
    args.k = 10;
    while ((op = getopt(argc, argv, "c:q:k:vmos:i")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'q':
                args.query_file = optarg;
                break;
            case 'k':
                args.k = std::strtoul(optarg, NULL, 10);
                break;
            case 'v':
                args.verbose = true;
                break;
            case 'm':
                args.multi_occ = true;
                break;
            case 'o':
                args.match_only = true;
                break;
            case 's':
                args.snippet_size = std::strtoul(optarg, NULL, 10);
                break;
            case 'i':
                args.intersection = true;
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
    if (args.collection_dir == "" || args.query_file == "") {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}



using idx_type = INDEX_TYPE;

const size_t buf_size = 1024 * 128;
char   buffer[buf_size];

template<typename X>
struct myline {
    static string parse(const string& str) {
        return str;
    }

    static vector<string> parse_multi(const string& str) {
        vector<string> terms;
        string cur_term;
        for (char c : str) {
            if (c == '\1') {
                terms.push_back(cur_term);
                cur_term.clear();
            } else {
                cur_term += c;
            }
        }
        terms.push_back(cur_term);
        return terms;
    }
};

template<>
struct myline<sdsl::int_alphabet_tag> {
    static vector<uint64_t> parse(const string& str) {
        vector<uint64_t> res;
        stringstream ss(str);
        uint64_t x;
        while (ss >> x) {
            res.push_back(x);
        }
        return res;
    }

    static vector<vector<uint64_t>> parse_multi(const string& str) {
        vector<vector<uint64_t>> terms;
        vector<uint64_t> cur_term;

        stringstream ss(str);
        uint64_t x;
        while (ss >> x) {
            if (x == 1) {
                terms.push_back(cur_term);
                cur_term.clear();
            } else {
                cur_term.push_back(x);
            }
        }
        terms.push_back(cur_term);
        return terms;
    }
};


template<typename T>
uint64_t get_input_size(string dir) {
    sdsl::int_vector<8> text;
    sdsl::load_from_file(text, dir + "/text_SURF.sdsl");
    return size_in_bytes(text);
}
template<>
uint64_t get_input_size<sdsl::int_alphabet_tag>(string dir) {
    sdsl::int_vector<> text;
    sdsl::load_from_file(text, dir + "/text_int_SURF.sdsl");
    return size_in_bytes(text);
}

int main(int argc, char* argv[])
{
    stick_this_thread_to_core(1);

    cmdargs_t args = parse_args(argc, argv);
    idx_type idx;

    using timer = chrono::high_resolution_clock;

    if (!args.verbose) {
        cout << "# collection_file = " << args.collection_dir << endl;
        cout << "# index_name = " << IDXNAME << endl;
    }
    auto cc = surf::parse_collection<idx_type::alphabet_category>(args.collection_dir);
    idx.load(cc);
    idx_type::topk_interface* topk = &idx;

    if (!args.verbose) {
        cout << "# pattern_file = " << args.query_file << endl;
        cout << "# doc_cnt = " << idx.doc_cnt() << endl;
        cout << "# word_cnt = " << idx.word_cnt() << endl;
        cout << "# k = " << args.k << endl;
        cout << "# match_only = " << args.match_only << std::endl;
        cout << "# multi_occ = " << args.multi_occ << std::endl;
    }
    ifstream in(args.query_file);
    if (!in) {
        cerr << "Could not load pattern file" << endl;
        return 1;
    }
    vector<string> queries;
    {   // Load all queries.
        string str;
        while (std::getline(in, str))
            queries.push_back(str);
    }

    using timer = chrono::high_resolution_clock;
    auto start = timer::now();

    vector<uint64_t> timings;  // in microseconds

    size_t sum = 0;
    size_t sum_fdt = 0;
    size_t sum_chars_extracted = 0;
    size_t q_len = 0;
    size_t q_cnt = 0;
    for (int run = 0; run < 1; ++run) {
        sum = 0;
        sum_fdt = 0;
        sum_chars_extracted = 0;
        q_len = 0;
        q_cnt = 0;
        for (size_t i = 0; i < queries.size(); ++i) {
            std::unique_ptr<idx_type::topk_interface::iter> res_it;
            idx_type::topk_interface::intersect_query intersect_query;
            if (args.intersection) {
                auto terms = myline<idx_type::alphabet_category>::parse_multi(
                        queries[i].c_str());
                for (const auto& term : terms) {
                    intersect_query.emplace_back(
                            term.data(), term.data() + term.size());
                    q_len += term.size();
                }

                ++q_cnt;
                start = timer::now();
                res_it = topk->topk_intersect(args.k, intersect_query,
                                              args.multi_occ, args.match_only);
            } else {
                auto query = myline<idx_type::alphabet_category>::parse(queries[i].c_str());
                q_len += query.size();
                ++q_cnt;
                start = timer::now();
                res_it = topk->topk(args.k, query.data(), query.data() + query.size(),
                                    args.multi_occ, args.match_only);
            }
            size_t x = 0;
            while (x < args.k && !res_it->done()) {
                ++x;
                sum_fdt += res_it->get().second;
                if (args.verbose) {
                    cout << q_cnt << ";" << x << ";" << res_it->get().first
                        << ";" << res_it->get().second << endl;
                }
                if (args.snippet_size != 0) {
                    auto snippet = res_it->extract_snippet(args.snippet_size);
                    sum_chars_extracted += snippet.size();
                    if (args.verbose) {
                        for (const auto c : snippet)
                            cout << c;
                        cout << endl;
                    }
                }
                if (x < args.k)
                    res_it->next();
            }
            sum += x;
            uint64_t msecs = chrono::duration_cast<chrono::microseconds>(
                    timer::now() - start).count();
            timings.push_back(msecs);
        }
    }


    if (!args.verbose) {
        sort(timings.begin(), timings.end());
        uint64_t qtime_median = timings[timings.size() / 2];
        uint64_t qtime_sum = accumulate(timings.begin(), timings.end(), 0);
        double qtime_avg = 1. * qtime_sum / timings.size();
        uint64_t qtime_max = *max_element(timings.begin(), timings.end());
        double qtime_sigma = 0;
        for (uint64_t t : timings)
            qtime_sigma += pow(t - qtime_avg, 2);
        qtime_sigma = sqrt(qtime_sigma / timings.size());

        cout << setprecision(5) << fixed;
        cout << "# query_len = " << q_len / q_cnt << endl;
        cout << "# queries = " << q_cnt << endl;
        cout << "# time_per_query_avg = " << qtime_avg << endl;
        cout << "# time_per_query_median = " << qtime_median << endl;
        cout << "# time_per_query_max = " << qtime_max << endl;
        cout << "# time_per_query_sigma = " << qtime_sigma << endl;
        auto doc_time = sum == 0 ? 0.0 : ((double)qtime_sum) / (sum * q_cnt);
        cout << "# time_per_doc = " << doc_time << endl;
        cout << "# check_sum = " << sum << endl;
        cout << "# check_sum_fdt = " << sum_fdt << endl;
        cout << "# index_size =  " << size_in_bytes(idx) << endl;
        cout << "# input_size = " <<
             get_input_size<idx_type::alphabet_category>(args.collection_dir) << endl;
        cout << "# sum_chars_extracted = " << sum_chars_extracted << endl;
    }
}
