
INCLUDE_FLAGS=-I/usr/include/eigen3 -I/usr/include/suitesparse

LIBS=-lg2o_core -lg2o_types_slam2d -lg2o_csparse_extension -lg2o_solver_csparse  -lg2o_stuff -lcxsparse -lg2o_opengl_helper

g2o_example: Uncertainty.o
	gcc Uncertainty.o -lstdc++ -lm ${LIBS} -o g2o_example
Uncertainty.o: Uncertainty.cpp
	gcc -c Uncertainty.cpp ${INCLUDE_FLAGS}

clean:
	rm Uncertainty.o g2o_example
