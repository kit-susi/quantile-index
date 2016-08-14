
#include "sdsl/config.hpp"
#include "sdsl/suffix_trees.hpp"
//#include "surf/indexes.hpp"
#include "surf/util.hpp"

using namespace sdsl;
using namespace std;

typedef struct cmdargs {
    std::string collection_dir;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -m\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    while ((op=getopt(argc,argv,"c:m:b")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
    if (args.collection_dir=="") {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

int main(int argc,char* const argv[])
{
    using clock = std::chrono::high_resolution_clock;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    sdsl::cache_config cc = surf::parse_collection<sdsl::byte_alphabet_tag>(args.collection_dir);
    // Input. 
    bit_vector H;
    int_vector<> P;
    int_vector<> weights;
    int_vector<> D;
    int_vector<> dup;
    { // Load files.
	    load_from_cache(H, surf::KEY_H, cc);
	    load_from_cache(P, surf::KEY_P, cc);
	    load_from_cache(weights, surf::KEY_WEIGHTS, cc);
	    load_from_cache(D, surf::KEY_DARRAY, cc);
	    load_from_cache(dup, surf::KEY_DUP, cc);
	    cout << "DUPSIZE: " << dup.size() << endl;
    }
    int_vector<> P_g(P);
    int_vector<> weights_g(weights);
    int_vector<> dup_g(dup);
    uint64_t output_pos = 0;
    select_support_mcl<> select_h(&H);
    for (uint64_t i = 2; i <= D.size(); ++i) {
	uint64_t ss = select_h(i-1);
	uint64_t ee  = select_h(i);
	if (ss < ee-1) { // Ignore empty intervals.
		uint64_t s = ss + 2 - i;
		uint64_t e = ee + 1 - i;
		// [s,e) interval.
		// Build set.
		std::unordered_map<uint64_t, uint64_t> dup_set;
		for (size_t j = s; j < e; ++j) {
			dup_set[dup[j]] = j;
		}
		uint64_t sa_pos = i;
		// Make order.
		while (!dup_set.empty()) {
			auto it = dup_set.find(D[sa_pos]);
			if (it != dup_set.end()) {
				// Insert found element next.
				P_g[output_pos] = P[it->second];
				weights_g[output_pos] = weights[it->second];
				dup_g[output_pos] = dup[it->second];
    				if (output_pos < 10) cout << dup_g[output_pos] << endl;
				output_pos++;
				dup_set.erase(it);
			}	
			if (sa_pos >= D.size())
				cout << "ERROR: sa_pos is out of bounds." << endl;
			++sa_pos;
		}
	}
    }
    store_to_cache(P_g, surf::KEY_P_G, cc);
    store_to_cache(weights_g, surf::KEY_WEIGHTS_G, cc);
    store_to_cache(dup_g, surf::KEY_DUP_G, cc);
    return EXIT_SUCCESS;
}
