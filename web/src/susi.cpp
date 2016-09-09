#include <Python.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <sys/types.h>
#include <dirent.h>

#include "failure.h"
#include "sdsl/sd_vector.hpp"

#define NAME IDX_NN
#define CSA_TYPE sdsl::csa_wt<sdsl::wt_huff<sdsl::hyb_vector<>>,32,32, sdsl::text_order_sa_sampling<>, sdsl::text_order_isa_sampling_support<>>
#define DF_TYPE surf::df_sada<sdsl::rrr_vector<63>,sdsl::rrr_vector<63>::select_1_type, sdsl::byte_alphabet_tag>
#define WTD_TYPE sdsl::wt_int<sdsl::bit_vector, sdsl::rank_support_v5<>, sdsl::select_support_scan<1>, sdsl::select_support_scan<0>>
#define KTWOTREAP_TYPE sdsl::k2_treap<2,sdsl::rrr_vector<63>>
#define INDEX_TYPE surf::idx_nn<CSA_TYPE,char,KTWOTREAP_TYPE>

// do we really need 3 levels here? https://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define XSTR(s) STR(s)
#define STR(s) #s
#define NAME_STR XSTR(NAME)

#include "surf/config.hpp"
#include "surf/idx_nn.hpp"

using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

namespace susi {
namespace python {
namespace {

using SusiIndex = INDEX_TYPE;

void free_index(PyObject* index_) {
    auto* index = reinterpret_cast<SusiIndex*>(PyCapsule_GetPointer(
            index_, "susi_index"));
    delete index;
}

#define PyBytesObject_SIZE (offsetof(PyBytesObject, ob_sval) + 1)

}  // namespace

static PyObject *
load(PyObject *self, PyObject *args) {
    const char *collection_dir_c;
    if (!PyArg_ParseTuple(args, "s", &collection_dir_c))
        return NULL;

    SusiIndex* idx = new SusiIndex();
    auto cc = surf::parse_collection<SusiIndex::alphabet_category>(collection_dir_c);
    idx->load(cc);

    return PyCapsule_New(reinterpret_cast<void*>(idx),
                         "susi_index", free_index);
}

static PyObject *
search(PyObject *self, PyObject *args) {
    PyObject* index_;
    const char *pattern;
    int pattern_size;
    int k;
    int snip_size;
    if (!PyArg_ParseTuple(args, "Os#ii",
                &index_, &pattern, &pattern_size, &k, &snip_size))
        return NULL;

    auto* index = reinterpret_cast<SusiIndex*>(PyCapsule_GetPointer(
            index_, "susi_index"));
    auto res = index->topk(k, pattern, pattern + pattern_size, false, false);

    auto* py_results = PyList_New(0);
    while (!res->done() && k--) {
        auto snippet = res->extract_snippet(snip_size);
        PyList_Append(py_results,
                PyTuple_Pack(3,
                    PyInt_FromLong(res->get().first),
                    PyFloat_FromDouble(res->get().second),
                    PyBytes_FromStringAndSize(
                        snippet.data(), snippet.size())));
        res->next();
    }
    return py_results;
}

static PyMethodDef Methods[] = {
    {"load",  load, METH_VARARGS,
    "Load a susi index from the given index directory."},
    {"search",  search, METH_VARARGS,
    "Use the given index to search for a pattern."},

    {NULL, NULL, 0, NULL} // sentinel
};

}  // namespace python
}  // namespace susi

PyMODINIT_FUNC
initsusi_native(void) {
    (void) Py_InitModule("susi_native", susi::python::Methods);
}
