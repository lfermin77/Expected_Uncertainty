
INCLUDE_FLAGS=-I/usr/include/eigen3 -I/usr/include/suitesparse

LIBS=-lg2o_core -lg2o_types_slam2d -lg2o_csparse_extension -lg2o_solver_csparse  -lg2o_stuff -lcxsparse -lg2o_opengl_helper

g2o_example: Expected_Several_Trajectories.o
	gcc Expected_Several_Trajectories.o -lstdc++ -g -lm ${LIBS} -o g2o_example
Expected_Several_Trajectories.o: Expected_Several_Trajectories.cpp
	gcc -c Expected_Several_Trajectories.cpp ${INCLUDE_FLAGS}

clean:
	rm Expected_Several_Trajectories.o g2o_example
