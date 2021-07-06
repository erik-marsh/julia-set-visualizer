CXX = g++
CXX_FLAGS = -std=c++17 -Wall -g

PROGNAME = julia
OUTFILE = out.tga
IMAGE_VIEWER = feh

.PHONY: clean all test profile

all: $(PROGNAME)

clean:
	rm -f $(OUTFILE)
	rm -f $(PROGNAME)
	rm -f profile*
	rm -f gmon.out

test: clean $(PROGNAME)
	./$(PROGNAME) 2048
	$(IMAGE_VIEWER) $(OUTFILE)

profile: main.cpp clean
	$(CXX) $(CXX_FLAGS) -g -pg $< -o $@
	./profile 4096
	gprof ./profile > profile.out

$(PROGNAME): main.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@
	