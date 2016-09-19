#include "sdsl/config.hpp"
#include "surf/util.hpp"
#include "sdsl/int_vector_buffer.hpp"
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
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -m\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -m <pattern length>        : the  pattern length.\n");
    fprintf(stdout,"  -x <number of patterns>    : generate x distinct patterns.\n");
    fprintf(stdout,"  -o <occurences>            : only return ngrams that are sampled at least once each X samples. 0 = disable\n");
    fprintf(stdout,"  -n <ngram samples>         : how many ngrams to sample if -p is given\n");
    fprintf(stdout,"  -s <random seed>\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.pat_len = 0;
    args.pat_cnt = 0;
    args.ngram_samples = 1<<13;
    args.min_sampling = 0;

    auto now = std::chrono::high_resolution_clock::now();
    srand(chrono::duration_cast<chrono::microseconds>(now.time_since_epoch()).count());
    args.seed = rand();

    while ((op=getopt(argc,argv,"c:m:x:o:s:n:")) != -1) {
        switch (op) {
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
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    if (args.collection_dir=="" || args.pat_len==0 || args.pat_cnt==0) {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

template <typename T, typename R>
void get_ngram(T& buf, R& dice, size_t n, string& res) {
    for(;;) {
        auto x = dice();
        bool valid = true;
        for (size_t j=x; j < x + n and valid; ++j){
            res[j-x] = buf[j];
            if ( buf[j] == 0 or buf[j] == 1 or buf[j] == '\n' or buf[j] == '\r')
                valid = false;
        }
        if (valid) return;
    }
}

int main(int argc,char* const argv[])
{
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);
//    std::cout<<"collection dir="<<args.collection_dir<<std::endl;
    std::string text_file = args.collection_dir+"/"+surf::TEXT_FILENAME_BYTE;
//    std::cout<<"> "<<text_file<<std::endl;
    sdsl::int_vector_buffer<8> text_buf(text_file, std::ios::in);
//    std::cout<<"n="<<text_buf.size()<<std::endl;

    if ( text_buf.size() < args.pat_len ){
        std::cout<<"ERROR: text.size()="<<text_buf.size()<<" < m=" <<args.pat_len << std::endl;
        return 1;
    }

    //cout << "seed = " << args.seed << endl;
    std::mt19937_64 rng(args.seed);
    std::uniform_int_distribution<uint64_t> distribution(0, text_buf.size()-args.pat_len);
    auto dice = bind(distribution, rng);

    string ngram(args.pat_len, '0');

    if (args.min_sampling) {
        unordered_map<string, size_t> ngram_counts(args.ngram_samples << 1);
        for (int i = args.ngram_samples; i--;) {
            get_ngram(text_buf, dice, args.pat_len, ngram);
            ngram_counts[ngram]++;
        }

        vector<string> ngrams;
        ngrams.reserve(args.pat_cnt);
        auto threshold = args.ngram_samples / args.min_sampling;
        for (const auto& it : ngram_counts)
            if (it.second >= threshold)
                ngrams.push_back(it.first);
        std::shuffle(ngrams.begin(), ngrams.end(), rng);

        if (args.pat_cnt > ngrams.size()) {
            std::cerr << "Not enough ngrams found" << endl;
            exit(1);
        }
        for (size_t i = 0; i < args.pat_cnt; ++i)
            cout << ngrams[i] << "\n";
    } else {
        unordered_set<string> results;
        while (results.size() < args.pat_cnt) {
            get_ngram(text_buf, dice, args.pat_len, ngram);
            results.insert(ngram);
        }
        for (const string& x: results)
            cout << x << "\n";
    }

    return EXIT_SUCCESS;
}


