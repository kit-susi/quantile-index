#include "sdsl/config.hpp"
#include "surf/util.hpp"
#include "sdsl/int_vector_buffer.hpp"
#include "sdsl/suffix_trees.hpp"
#include <algorithm>
#include <chrono>
#include <map>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace std;

typedef struct cmdargs {
    std::string collection_dir;
    size_t pat_len;
    size_t pat_cnt;
    size_t min_sampling;
    size_t seed;
    size_t ngram_samples;
    bool int_collection;
} cmdargs_t;

void
usage(char* program)
{
    printf("%s -c <collection directory> -m\n", program);
    printf("where\n");
    printf("  -c <collection directory>  : the directory the collection is stored.\n");
    printf("  -m <pattern length>        : the  pattern length.\n");
    printf("  -x <number of patterns>    : generate x distinct patterns.\n");
    printf("  -o <occurences>            : only return ngrams that are sampled at least once each X samples. 0 = disable (assertion o >= n)\n");
    printf("  -n <ngram samples>         : how many ngrams to sample if -o is given\n");
    printf("  -s <random seed>\n");
    exit(EXIT_FAILURE);
};

cmdargs_t
parse_args(int argc, char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.pat_len = 0;
    args.pat_cnt = 0;
    args.ngram_samples = 1 << 13;
    args.min_sampling = 0;
    args.int_collection = false;

    auto now = std::chrono::high_resolution_clock::now();
    srand(chrono::duration_cast<chrono::microseconds>(now.time_since_epoch()).count());
    args.seed = rand();

    while ((op = getopt(argc, argv, "ic:m:x:o:s:n:")) != -1) {
        switch (op) {
            case 'i':
                args.int_collection = true;
                break;
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'm':
                args.pat_len = std::stoull(std::string(optarg));
                break;
            case 'x':
                args.pat_cnt = std::stoull(std::string(optarg));
                break;
            case 'o':
                args.min_sampling = std::stoull(std::string(optarg));
                break;
            case 's':
                args.seed = std::stoull(std::string(optarg));
                break;
            case 'n':
                args.ngram_samples = std::stoull(std::string(optarg));
                break;
            case '?':
            default:
                usage(argv[0]);
        }
    }
    if (args.collection_dir == "" || args.pat_len == 0 || args.pat_cnt == 0) {
        std::cerr << "Missing command line parameters.\n";
        usage(argv[0]);
    } else if (args.ngram_samples < args.min_sampling) {
        std::cerr << "o should be greater then n.\n";
        usage(argv[0]);
    }
    return args;
}

template <typename T, typename R, typename S>
void get_ngram(T& buf, R& dice, size_t n, S& res) {
    for (;;) {
        auto x = dice();
        bool valid = true;
        for (size_t j = x; j < x + n and valid; ++j) {
            res[j - x] = buf[j];
            if (buf[j] == 0 or buf[j] == 1 or buf[j] == '\n' or buf[j] == '\r')
                valid = false;
        }
        if (valid) return;
    }
}

class string_vector_hash_fn {
public:
    std::size_t operator()(std::vector<int> const& vec) const {
        std::size_t ret = 0;
        for (auto & i : vec) {
            ret ^= std::hash<uint32_t>()(i);
        }
        return ret;
    }

    std::size_t operator()(std::string const& s) const {
        std::hash<string> h;
        return h(s);
    }
};


template<typename S, typename T>
int find_ngrams(T& text_buf, cmdargs_t& args) {
    if (text_buf.size() < args.pat_len) {
        std::cout << "ERROR: text.size()=" << text_buf.size() << " < m=" << args.pat_len << std::endl;
        return 1;
    }

    //cout << "seed = " << args.seed << endl;
    std::mt19937_64 rng(args.seed);
    std::uniform_int_distribution<uint64_t> distribution(0, text_buf.size() - args.pat_len);
    auto dice = bind(distribution, rng);

    S ngram(args.pat_len, '0');

    if (args.min_sampling) {
        string_vector_hash_fn h;
        unordered_map<size_t, S> ngram_storage(args.ngram_samples << 1);
        unordered_map<size_t, size_t> ngram_counts(args.ngram_samples << 1);
        for (int i = args.ngram_samples; i--;) {
            get_ngram(text_buf, dice, args.pat_len, ngram);
            ngram_counts[h(ngram)]++;
            ngram_storage[h(ngram)] = ngram;
        }

        vector<S> ngrams;
        ngrams.reserve(args.pat_cnt);
        auto threshold = args.ngram_samples / args.min_sampling;
        for (const auto & it : ngram_counts)
            if (it.second >= threshold)
                ngrams.push_back(ngram_storage[it.first]);
        std::shuffle(ngrams.begin(), ngrams.end(), rng);

        if (args.pat_cnt > ngrams.size()) {
            std::cerr << "Not enough ngrams found" << endl;
            exit(1);
        }
        for (size_t i = 0; i < args.pat_cnt; ++i)
            for (size_t j = 0; j < ngrams[i].size(); ++j)
                cout << ngrams[i][j] << ((j + 1 == ngrams[i].size()) ? "\n" : (args.int_collection ? " " : ""));
    } else {
        string_vector_hash_fn h; vector<S> results;
        unordered_map<size_t, bool> vis;
        while (results.size() < args.pat_cnt) {
            get_ngram(text_buf, dice, args.pat_len, ngram);
            if (!vis[h(ngram)]) {
                results.push_back(ngram); vis[h(ngram)] = true;
            }
        }
        for (const S & x : results)
            for (size_t j = 0; j < x.size(); ++j)
                cout << x[j] << ((j + 1 == x.size()) ? "\n" : (args.int_collection ? " " : ""));
    }

    return EXIT_SUCCESS;
}

int main(int argc, char* const argv[])
{
    /* parse command line */
    cmdargs_t args = parse_args(argc, argv);
    if (args.int_collection) {
        std::string text_file = args.collection_dir + "/" + surf::TEXT_FILENAME;
        sdsl::int_vector_buffer<> text_buf(text_file, std::ios::in);
        return find_ngrams<vector<int>>(text_buf, args);
    }
    else {
        std::string text_file = args.collection_dir + "/" + surf::TEXT_FILENAME_BYTE;
        sdsl::int_vector_buffer<8> text_buf(text_file, std::ios::in);
        return find_ngrams<string>(text_buf, args);
    }

}


