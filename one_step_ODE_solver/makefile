CXXFLAGS += -std=c++17 -pedantic -pedantic-errors -Wall
INCL_DIR = include
SRC_DIR = src

main: main.o jacobian.o SLAE_solver.o SNAE_solver.o operations.o one_step_ODE_solver.o
	$(CXX) -o main main.o jacobian.o SLAE_solver.o SNAE_solver.o operations.o one_step_ODE_solver.o
	
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -I$(INCL_DIR) -c main.cpp
	
jacobian.o: $(SRC_DIR)/jacobian.cpp
	$(CXX) $(CXXFLAGS) -I$(INCL_DIR) -c $(SRC_DIR)/jacobian.cpp
	
operations.o: $(SRC_DIR)/operations.cpp
	$(CXX) $(CXXFLAGS) -I$(INCL_DIR) -c $(SRC_DIR)/operations.cpp
	
SNAE_solver.o: $(SRC_DIR)/SNAE_solver.cpp
	$(CXX) $(CXXFLAGS) -I$(INCL_DIR) -c $(SRC_DIR)/SNAE_solver.cpp
	
SLAE_solver.o: $(SRC_DIR)/SLAE_solver.cpp
	$(CXX) $(CXXFLAGS) -I$(INCL_DIR) -c $(SRC_DIR)/SLAE_solver.cpp
	
one_step_ODE_solver.o: $(SRC_DIR)/one_step_ODE_solver.cpp
	$(CXX) $(CXXFLAGS) -I$(INCL_DIR) -c $(SRC_DIR)/one_step_ODE_solver.cpp 
	
clean:
	rm -rf *.o main
