#include <sdsl/int_vector.hpp>
#include <iostream>
#include <sstream>
#include <unordered_set>

using namespace sdsl;
using namespace std;

bool int_collection = false;
size_t num_docs = 0;
size_t doc_len = 0;
string infile;
string outdir = ".";

template <typename T>
void process(T& text) {
    uint64_t num_docs = 0;
    unordered_set<uint64_t> in_use;
    for (size_t i = 0; i < text.size(); ++i) {
        if (text[i] == 1)
            num_docs++;
        in_use.insert(text[i]);
    }
    cout << text.size() << " " << num_docs << " " << text.size() / num_docs << " "
         << in_use.size() << " " << text.size() / 1024/1024<< endl;
}

void usage(const string& prog) {
    cout << "Usage: " << prog << " [-i] file_name" << endl;
    exit(EXIT_FAILURE);
}

void parse_opts(int argc, char* const argv[]) {
    int op;
    while ((op = getopt(argc, argv, "i:")) != -1) {
        switch (op) {
            case 'i':
                int_collection = true;
                break;
            case '?':
            default:
                usage(argv[0]);
        }
    }
    if (argc - optind != 1) usage(argv[0]);
    infile = argv[optind];
}

int main(int argc, char* argv[]) {
    parse_opts(argc, argv);
    if (int_collection) {
        int_vector<> text;
        load_from_file(text, infile);
        process(text);
    } else {
        int_vector<8> text;
        load_vector_from_file(text, infile, 1);
        process(text);
    }
}
