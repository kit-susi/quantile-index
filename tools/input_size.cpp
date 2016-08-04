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
    load_from_file(text, argv[1]);
    cout << size_in_bytes(text) << endl;
    return 0;
}
