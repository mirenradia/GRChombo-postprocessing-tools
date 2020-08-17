# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

GRCHOMBO_SOURCE ?= ../../../GRChombo/Source

ebase := Norm

LibNames := BaseTools

src_dirs := ../Source \
            ../Source/TimeData \
            $(GRCHOMBO_SOURCE)/AMRInterpolator \
            $(GRCHOMBO_SOURCE)/utils

# Filesystem library is a C++17 feature
XTRACXXFLAGS := $(XTRACXXFLAGS) -std=c++17
#export XTRACXXFLAGS

# Older versions of GCC require explicit linking with the Filesystem library
XTRALDFLAGS := $(XTRALDFLAGS) -lstdc++fs
#export XTRALDFLAGS

include $(CHOMBO_HOME)/mk/Make.test
