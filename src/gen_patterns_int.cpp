#include "sdsl/config.hpp"
#include "surf/util.hpp"
#include "sdsl/int_vector_buffer.hpp"
#include <chrono>
#include <random>
using namespace std;

//const size_t buf_size = 1024*1024;
//char buffer[buf_size];

typedef struct cmdargs {
    std::string collection_dir;
    size_t pat_len;
    size_t pat_cnt;
    size_t seed;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -m\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -m <pattern length>        : the  pattern length.\n");
    fprintf(stdout,"  -x <number of patterns>    : generate x patterns.\n");
    fprintf(stdout,"  -s <random seed>\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.pat_len  = 0;
    args.pat_cnt  = 0;

    auto now = std::chrono::high_resolution_clock::now();
    srand(chrono::duration_cast<chrono::microseconds>(now.time_since_epoch()).count());
    args.seed = rand();

    while ((op=getopt(argc,argv,"c:m:x:")) != -1) {
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
            case 's':
                args.seed = std::stoull(std::string(optarg));
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
    if (args.collection_dir=="" || args.pat_len==0 || args.pat_cnt==0) {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

int main(int argc,char* const argv[])
{
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);
//    std::cout<<"collection dir="<<args.collection_dir<<std::endl;
    std::string text_file = args.collection_dir+"/"+surf::TEXT_FILENAME;
//    std::string dict_file = args.collection_dir+"/"+surf::DICT_FILENAME;
//    std::cout<<"> "<<text_file<<std::endl;
    sdsl::int_vector_buffer<> text_buf(text_file, std::ios::in);
//    std::cout<<"n="<<text_buf.size()<<std::endl;

    if ( text_buf.size() < args.pat_len ){
        std::cout<<"ERROR: text.size()="<<text_buf.size()<<" < m=" <<args.pat_len << std::endl;
        return 1;
    }

    std::mt19937_64 rng(args.seed);
    std::uniform_int_distribution<uint64_t> distribution(0, text_buf.size()-args.pat_len);
    auto dice = bind(distribution, rng);

    for (size_t i=0; i<args.pat_cnt; ++i){ 
        auto x = dice();
        bool valid = true;
        for (size_t j=x; j < x + args.pat_len and valid; ++j){
            if ( text_buf[j] == 0 or text_buf[j] == 1)
                valid = false;
        }
        if ( valid ){
            for (size_t j=x; j < x + args.pat_len; ++j){
                std::cout<<text_buf[j]<<" ";
            }
            std::cout<<std::endl;
        } else {
            --i;
        }
    }
    return EXIT_SUCCESS;
}
