CXX = g++
CXX_FLAGS = -std=c++17 -Wall -lpthread

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
	./$(PROGNAME) 2048 -0.8 0.156 1
	$(IMAGE_VIEWER) $(OUTFILE)

profile: main.cpp clean
	$(CXX) $(CXX_FLAGS) -g -pg $< -o $@
	./profile 4096 -0.8 0.156 1
	gprof ./profile > profile.out

profile-long: main.cpp clean
	$(CXX) $(CXX_FLAGS) -g -pg $< -o $@
	./profile 4096 -0.8 0.156 1
	gprof ./profile > profile-1t.out
	./profile 4096 -0.8 0.156 2
	gprof ./profile > profile-2t.out
	./profile 4096 -0.8 0.156 3
	gprof ./profile > profile-3t.out
	./profile 4096 -0.8 0.156 4
	gprof ./profile > profile-4t.out

$(PROGNAME): main.cpp
	$(CXX) $(CXX_FLAGS) $< -o $@
	