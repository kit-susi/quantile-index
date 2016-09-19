#include <sdsl/int_vector.hpp>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <map>

using namespace boost;
using namespace sdsl;
using namespace std;

const string SEP = "\1";


bool valid(const string& str) {
    if (str.size() < 2)
        return false;
    return true;
}

string space() {
    string res = "";
    for (int a = 0; a < 256; ++a) {
        char c(a);
        if (std::isspace(c))
            res += c;
    }
    res += '\1';
    return res;
}

string punct() {
    string res = "";
    for (int a = 0; a < 256; ++a) {
        char c(a);
        if (std::ispunct(c))
            res += c;
    }
    return res;
}

string outdir = ".";
string infile;
bool input_is_sdsl = false;

void usage(string prog) {
    cout << "Usage: " << prog << " [-o outdir] [-s] file_name" << endl;
    cout << "  Takes a text file containing documents seperated with 1 (binary)." << endl;
    cout << "  Stores parsed collection in `text_int_SURF.sdsl` and `dict.txt`." << endl;
    cout << "  If -s is given, input is assumed to be a text surf collection " << endl;
    cout << "  (i.e. a text_SURF.sdsl file)" << endl;
    exit(EXIT_FAILURE);
}

void parse_opts(int argc, char* argv[]) {
    int op;
    while ((op = getopt(argc, argv, "o:s")) != -1) {
        switch (op) {
            case 'o':
                outdir = optarg;
                break;
            case 's':
                input_is_sdsl = true;
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

    //int_vector<> text;
    //load_from_file(text, argv[1]);
    //text.resize(text.size()+1);
    //text[text.size()-1]=0;
    //store_to_file(text, "text_int_SURF.sdsl");
    typedef boost::tokenizer<boost::char_separator<char>,
    std::istreambuf_iterator<char> > tokenizer;
    boost::char_separator<char> sep(punct().c_str(), space().c_str());
    unordered_map<string, uint64_t> dictionary;
    uint64_t num_docs = 0;
    uint64_t total_size = 0;
    {   // Build dictionary.
        ifstream ifile(infile);
        if (input_is_sdsl) {
            // skip int vector header
            char throwaway[8];
            ifile.read(throwaway, 8);
        }
        istreambuf_iterator<char> file_iter(ifile);
        istreambuf_iterator<char> end_of_stream;
        tokenizer tokens(file_iter, end_of_stream, sep);
        for (const auto& word : tokens) {
            if (word == SEP) {
                num_docs++;
                total_size++;
            } else if (valid(word)) {
                if (dictionary.count(word) == 0) {
                    uint64_t id = dictionary.size();
                    dictionary[word] = id;
                }
                total_size++;
            }
        }
    }
    cout << "num_docs " << num_docs << endl;
    cout << "dict_size " << dictionary.size() << endl;
    cout << "total_size " << total_size << endl;
    int_vector<> text;

    text.resize(total_size + 1);
    text[text.size() - 1] = 0;

    vector<uint64_t> word_count(dictionary.size(), 0);
    {   // Count words and built integer sequence.
        ifstream ifile(infile);
        if (input_is_sdsl) {
            // skip int vector header
            char throwaway[8];
            ifile.read(throwaway, 8);
        }
        istreambuf_iterator<char> file_iter(ifile);
        istreambuf_iterator<char> end_of_stream;
        tokenizer tokens(file_iter, end_of_stream, sep);
        uint64_t i  = 0;
        for (const auto& word : tokens) {
            if (word == SEP) {
                text[i] = 1;
                ++i;
            } else if (valid(word)) {
                if (dictionary.count(word)) {
                    text[i] = dictionary[word] + 2; // +2 because \0 \1 seperators.
                    ++word_count[dictionary[word]];
                }
                ++i;
            }
        }
    }
    for (uint64_t x: text)
        cout << x << endl;

    util::bit_compress(text);
    //for (const auto& w : text) cout << w << " "; cout << endl;
    store_to_file(text, outdir + "/text_int_SURF.sdsl");
    // Sort dict.
    vector<pair<string, int>> dict(dictionary.begin(), dictionary.end());
    sort(dict.begin(), dict.end(),  [&](const auto & a, const auto & b) {
        return word_count[a.second] > word_count[b.second];
    });
    {   // Output dictionary.
        fstream fout(outdir + "/dict.txt", ios::out);
        for (const auto& entry : dict) {
            fout << entry.first << " " << entry.second + 2 << " " <<
                 word_count[entry.second] << endl;
        }
    }
}
