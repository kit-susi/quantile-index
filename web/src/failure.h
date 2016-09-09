#pragma once

#ifdef USE_PYTHON
#  include <Python.h>
#else
#  define PyErr_SetString(a, b) do {} while(0)
#endif
#include <sstream>

#define PY_ERROR(ty, msg) \
    do { \
        std::stringstream ss; ss << msg; \
        PyErr_SetString(ty, ss.str().c_str()); \
        return false; \
    } while(0)

#define PY_ASSERT(expr, msg) \
    do { \
        if (!(expr)) { \
            PY_ERROR(PyExc_AssertionError, #expr << ": " << msg); \
        } \
    } while (0)

#define PROPAGATE_FAILURE(x) if (!(x)) return false
#define RETURN_SUCCESS return true
