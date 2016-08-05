#include <sdsl/int_vector.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;


int main(int argc, char* argv[]){
    if ( argc < 2 ){
        cout<<"Usage: "<<argv[0]<<" file_name"<<endl;
        cout<<" Prints the size in bytes of the passed int vector file."<<endl;
    }
    int_vector<> text;
    //int_vector<8> text;
    load_from_file(text, argv[1]);
    cout << size_in_bytes(text) << endl;
    uint64_t sigma = 0;
    int doc_count = 0;
    for (uint64_t i = 0; i < text.size(); ++i) {
	uint64_t n = text[i];
	sigma = std::max(sigma, n);
	if (n == 1) 
		doc_count++;
    }
    cout << "sigma: " << sigma << endl;
    cout << "doc_count: " << doc_count << endl;
    return 0;
}
