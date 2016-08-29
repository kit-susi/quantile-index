#include <sdsl/int_vector.hpp>
#include <iostream>
#include <sstream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
    if (argc < 2 || argc > 3){
        cout<<"Usage: "<<argv[0]<<" file_name [num_docs]"<<endl;
        cout<<" Takes byte character text and generates int_vector<8>"<<endl;
        cout<<" stored in `text_SURF.sdsl`. If num_docs is provided, "<<endl;
        cout<<" only the first num_docs documents are used."<<endl;
        return 1;
    }
    int_vector<8> text;
    load_vector_from_file(text, argv[1], 1);

    if (argc > 2) {
        int num_docs;
        stringstream(string(argv[2])) >> num_docs;
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
    store_to_file(text, "text_SURF.sdsl");
}
