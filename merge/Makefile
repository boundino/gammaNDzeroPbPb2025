CXX = g++
# CXXFLAGS += -O2 -Wall -Werror -Wextra
CXXFLAGS += -std=c++20
ROOTFLAGS := `root-config --cflags --libs` -lEG -I$(HOME)

BUILDDIR = ./build

SRCS  = $(wildcard *.cc)
EXES  = $(patsubst %.cc,%.exe,$(SRCS))
DEPS  = $(patsubst %.cc,$(BUILDDIR)/%.d,$(SRCS))

.PHONY: all clean

all: $(EXES)

%.exe: %.cc
	@mkdir -p $(BUILDDIR)/$(@D)
	$(CXX) $(ROOTFLAGS) $(CXXFLAGS) -MMD -MF $(BUILDDIR)/$(@D)/$(*F).d $< -o $@

clean:
	@$(RM) $(EXES) $(DEPS)
	@rm -rf $(BUILDDIR)/*

-include $(DEPS)
