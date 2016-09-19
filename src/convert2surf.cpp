#include <sdsl/int_vector.hpp>
#include <iostream>
#include <sstream>

using namespace sdsl;
using namespace std;

bool int_collection = false;
size_t num_docs = 0;
size_t doc_len = 0;
string infile;
string outdir = ".";

template <typename T>
void process(T& text) {
    if (doc_len) {
        T text_new(text.size());
        cout << "Using " << doc_len << " terms of first " << num_docs << " documents" << endl;
        size_t doc_start = 0;
        size_t new_offset = 0;
        while (num_docs > 0 && doc_start < text.size()) {
            // copy document from text to text_new, but only first doc_len
            // chars
            size_t i = 0;
            for (; i < doc_len && text[doc_start + i] != 1; ++i)
                text_new[new_offset++] = text[doc_start + i];
            text_new[new_offset++] = 1;

            doc_start += i;
            // find actual end of document
            while (doc_start < text.size() && text[doc_start] != 1)
                ++doc_start;
            ++doc_start;

            --num_docs;
        }
        text_new.resize(new_offset);
        std::swap(text_new, text);
    } else if (num_docs) {
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
    text.resize(text.size() + 1);
    text[text.size() - 1] = 0;
}

void usage(const string& prog) {
    cout << "Usage: " << prog << " [-i] [-d num_docs] [-l doc_len] [-o outdir] file_name" << endl;
    exit(EXIT_FAILURE);
}

void parse_opts(int argc, char* const argv[]) {
    int op;
    while ((op = getopt(argc, argv, "id:l:o:")) != -1) {
        switch (op) {
            case 'd':
                num_docs = stoull(string(optarg));
                break;
            case 'l':
                doc_len = stoull(string(optarg));
                break;
            case 'o':
                outdir = optarg;
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
    if (int_collection) {
        int_vector<> text;
        load_from_file(text, infile);
        process(text);
        store_to_file(text, outdir + "/text_int_SURF.sdsl");
    } else {
        int_vector<8> text;
        load_vector_from_file(text, infile, 1);
        process(text);
        store_to_file(text, outdir + "/text_SURF.sdsl");
    }
}
