#include <sdsl/int_vector.hpp>
#include <iostream>
#include <sstream>

using namespace sdsl;
using namespace std;

template <typename T>
void finalize(T& text, size_t num_docs=0) {
    if (num_docs) {
        cout << "Using first " << num_docs << " documents" << endl;
        size_t i = 0;
        while (num_docs > 0 && i < text.size()) {
            --num_docs;
            while (i < text.size() && text[i] != 1)
                ++i;
            ++i;
        }
        cout << "Truncating at index " << i << " / " << text.size() << endl;
        text.resize(i);
    }
    text.resize(text.size()+1);
    text[text.size()-1]=0;
}

void convert_text(const string& infile, size_t num_docs=0) {
    int_vector<8> text;
    load_vector_from_file(text, infile, 1);
    finalize(text);
    store_to_file(text, "text_SURF.sdsl");
}

void convert_int(const string& infile, int num_docs=0) {
    int_vector<> text;
    load_from_file(text, infile);
    finalize(text);
    store_to_file(text, "text_int_SURF.sdsl");
}

void usage(const string& prog) {
    cout<<"Usage: "<<prog<<" [-i] [-d num_docs] file_name"<<endl;
    exit(EXIT_FAILURE);
}

bool int_collection = false;
size_t num_docs = 0;
string infile;

void parse_opts(int argc, char* const argv[]) {
    int op;
    while ((op = getopt(argc, argv, "id:")) != -1) {
        switch (op) {
            case 'd':
                num_docs = stoull(string(optarg));
                break;
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
    if (int_collection)
        convert_int(infile, num_docs);
    else
        convert_text(infile, num_docs);
}
