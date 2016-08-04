#include <sdsl/int_vector.hpp>
#include <iostream>
#include<boost/tokenizer.hpp>
#include <map>

using namespace boost;
using namespace sdsl;
using namespace std;

const string SEP = "\1";

int main(int argc, char* argv[]){
    if ( argc < 2 ){
        cout<<"Usage: "<<argv[0]<<" file_name"<<endl;
        cout<<" Takes a text file containing documents seperated with 1 (binary)."<<endl;
        cout<<" Stores parsed collection in `text_int_SURF.sdsl` and `dict.txt`"<<endl;
    }
    //int_vector<> text;
    //load_from_file(text, argv[1]);
    //text.resize(text.size()+1);
    //text[text.size()-1]=0;
    //store_to_file(text, "text_int_SURF.sdsl");
    typedef boost::tokenizer<boost::char_separator<char>,
	    std::istreambuf_iterator<char> > tokenizer;
    unordered_map<string,uint64_t> dictionary; 
    uint64_t num_docs = 1;
    uint64_t total_size = 0;
    { // Build dictionary.
	    ifstream ifile(argv[1]);
	    istreambuf_iterator<char> file_iter(ifile);
	    istreambuf_iterator<char> end_of_stream;
	    tokenizer tokens(file_iter,end_of_stream);
	    for (const auto& word : tokens) {
		    if (word == SEP) {
			num_docs++;
		    } else if(dictionary.count(word) == 0) {
			uint64_t id = dictionary.size();
			dictionary[word] = id;
		    }
		    total_size++;
	    }
    }
    cout << "num_docs " << num_docs << endl;
    cout << "dict_size " << dictionary.size() << endl;
    cout << "total_size " << total_size << endl;
    int_vector<> text;
    text.resize(total_size+1);
    text[total_size]=0;
    vector<uint64_t> word_count(dictionary.size(), 0);
    { // Count words and built integer sequence.
    	    ifstream ifile(argv[1]);
	    istreambuf_iterator<char> file_iter(ifile);
	    istreambuf_iterator<char> end_of_stream;
	    tokenizer tokens(file_iter,end_of_stream);
	    uint64_t i  = 0;
	    for (const auto& word : tokens) {
		    if (word == SEP) {
		  	text[i] = 1;	  	
		    } else if(dictionary.count(word)) {
			text[i] = dictionary[word]+2; // +2 because \0 \1 seperators.
			++word_count[dictionary[word]];
		    }
		    ++i;
	    }
    }
    //for (const auto& w : text) cout << w << " "; cout << endl;
    store_to_file(text, "text_int_SURF.sdsl");
    { // Output dictionary.
	fstream fout("dict.txt", ios::out);
	for (const auto& entry : dictionary) {
		fout << entry.first << " " << entry.second + 2 << " " <<
		word_count[entry.second] << endl;	
	}
    }
}
