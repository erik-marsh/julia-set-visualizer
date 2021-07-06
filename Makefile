CXX = g++
CXX_FLAGS = -std=c++17 -Wall

.PHONY: clean all

all: julia

clean:
	rm -f out.tga
	rm -f julia

test-noimg: clean julia
	./julia 2048
	feh out.tga

test: test-noimg
	feh out.tga

julia: main.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@