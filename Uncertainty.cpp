/*
 * test.cpp
 *
 *  Created on: Sep 31, 2014
 *      Author: Leonardo
 */




#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/factory.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/core/sparse_block_matrix.h>

#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam2d/edge_se2.h>

#include <Eigen/Cholesky>

typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  		SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> 	SlamLinearSolver;




float compute_current_covariance(g2o::SparseOptimizer& optimizer, int Current_Node){

	g2o::OptimizableGraph::Vertex* v0;// = new g2o::VertexSE2;
	v0 = optimizer.vertex(Current_Node);
	g2o::SparseBlockMatrix<Eigen::MatrixXd> spinv;	

	optimizer.computeMarginals(spinv, v0);	
	int block_index = v0->colInHessian();


	std::cout << "asdfa  "<<std::endl;
			
	if ( block_index >= 0){
		Eigen::MatrixXd Remaped_Matrix(3,3);		
		for(int i=0 ; i < 3 ; i ++){
			for (int j=0 ; j < 3 ; j ++){
				Remaped_Matrix(i,j) = (*spinv.block(block_index/3,block_index/3))(i,j);
			}
		}
		return Remaped_Matrix.determinant();
	}
	else{
		return 0;
	}
	std::cout << "asdfa  "<<Remaped_Matrix<<std::endl;
}




int main(int argc, char** argv)
{

	g2o::EdgeSE2 edge;
	g2o::SparseOptimizer optimizer;

///////////// Initializing g2o related

	SlamLinearSolver* linearSolver = new SlamLinearSolver();
	linearSolver->setBlockOrdering(false);
	SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
	g2o::OptimizationAlgorithmGaussNewton* solverGauss   = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);

	optimizer.setAlgorithm(solverGauss);
	optimizer.load(argv[1]);
//	optimizer.vertex(2700)->setFixed(true);
	optimizer.vertex(0)->setFixed(true);
	optimizer.setVerbose(true);
	std::cerr<< "Antes de inicializar" <<std::endl;
//	optimizer.setVerbose(false);
	optimizer.initializeOptimization();
	std::cerr<< "Despues de inicializar" <<std::endl;
	
	optimizer.optimize(2);	
	
		std::cout << "asdfa  "<<std::endl;
		
	compute_current_covariance(optimizer, 100)

}

