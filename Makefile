CC = g++
CFLAGS = -g -I /usr/include/eigen3/ -I. -fopenmp -O3 -fdiagnostics-color=always -std=c++17

SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.d)
EXEC = main

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -L/usr/local/lib -lnlopt 

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	$(CC) -MM $< -MF $(@:.o=.d)

-include $(DEPS)

clean:
	rm -f $(OBJS) $(DEPS) $(EXEC)
