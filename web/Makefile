PY_SRCS := $(wildcard src/*.cpp)
PY_OBJS := $(subst src,build/src_python,$(patsubst %.cpp,%.o,$(PY_SRCS)))

SUSI := ..

DSS := ../external/sdsl-lite/external/libdivsufsort
DSS_BUILD := ../build/external/sdsl-lite/external/libdivsufsort
DSS_SRCS := $(wildcard $(DSS)/lib/*.c)
DSS_OBJS := $(subst $(DSS)/lib,build/dss,$(patsubst %.c,%.o,$(DSS_SRCS)))

SDSL := ../external/sdsl-lite
SDSL_BUILD := ../build/external/sdsl-lite
SDSL_SRCS := $(wildcard $(SDSL)/lib/*.cpp)
SDSL_OBJS := $(subst $(SDSL)/lib,build/sdsl,$(patsubst %.cpp,%.o,$(SDSL_SRCS)))


INCLUDES := \
    -I$(SUSI)/include \
    -I$(SDSL)/include \
    -I$(SDSL_BUILD)/include \
    -I$(DSS)/include \
    -I$(DSS_BUILD)/include \
    -Isrc \

CFLAGS := -Wall -g -O3 $(INCLUDES)

CPPFLAGS := $(CFLAGS) -std=c++14 -Wall -Wextra -Wpedantic -Wno-unused-parameter

PYFLAGS := \
    -I/usr/include/python2.7 \
    -DUSE_PYTHON \
    -fPIC \

.PHONY: all install clean pyobjs pybuild test

all: pybuild test

install: pybuild
	python2 setup.py install

test: pybuild
	PYTHONPATH=$(wildcard build/lib.*) python2 -m unittest susi.query_parser

pyobjs: $(PY_OBJS) $(SDSL_OBJS) $(DSS_OBJS)

pybuild: pyobjs
	rm -rf build/lib.*
	python2 setup.py build

build/src_python/%.o: src/%.cpp $(DSS_TARGET)
	mkdir -p build/src_python
	g++ -fPIC $(CPPFLAGS) $(PYFLAGS) -c $< -o $@

build/dss/%.o: $(DSS)/lib/%.c
	mkdir -p build/dss
	gcc -fPIC $(CFLAGS) -DHAVE_CONFIG_H=1 -c $< -o $@

build/sdsl/%.o: $(SDSL)/lib/%.cpp
	mkdir -p build/sdsl
	g++ -fPIC $(CPPFLAGS) -c $< -o $@

clean:
	rm -rf build
