#run:
#	source ~/srir/source_bash.sh && mpicxx -ggdb -std=c++1y map.cpp dijkstra.cpp main.cpp -DDEBUG -o main && mpiexec -n 2 ./main testcases/testcase0

# add -DDEBUG to turn on logging in all: line below (right before "-o main")
all:
	mpicxx -ggdb -std=c++1y map.cpp dijkstra.cpp main.cpp -o main
	
test0:
	mpiexec -n 3 ./main testcases/testcase0

test1:
	mpiexec -n 3 ./main testcases/testcase1

test3:
	mpiexec -n 4 ./main testcases/testcase3
