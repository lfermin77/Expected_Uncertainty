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


class Hypothesis{

	public:
	float probability_;
	Eigen::Matrix3d Relative_Uncertainty_;
	Eigen::Matrix3d Absolute_Uncertainty_;
	Eigen::Matrix3d First_Node_Uncertainty_;
	int last_loop_closure_;
	      Hypothesis(int a){
			  last_loop_closure_=0;
		  }
		 ~Hypothesis()
		{
//	do nothing
		}
	      
};

class Path_Quantization{

	public:

	Eigen::Matrix3d Last_Node_Uncertainty_to_First;
	float Last_Node_Uncertainty_to_First_Determinant;			

	Eigen::Matrix3d Expected_Uncertainty_First_Node;
	float Expected_Uncertainty_First_Node_Determinant;			


	Eigen::Matrix3d Expected_Uncertainty_Absolute;
	float Expected_Uncertainty_Absolute_Determinant;			

	Eigen::Matrix3d Last_node_Uncertainty;
	float Last_node_Uncertainty_Determinant;			

	float Last_node_probability;		
	
	float D_Optimal_Criteria;		
	float Work;
	
	float Distance;		
	
	int last_loop_closure_;
	      Path_Quantization(int a){
			  last_loop_closure_=0;
		  }
		 ~Path_Quantization()
		{
//	do nothing
		}
		
		void print_path_properties(){
				std::cout << "Last_Node_Uncertainty_to_First   "<< std::endl << Last_Node_Uncertainty_to_First_Determinant << std::endl;
				std::cout << "Expected_Uncertainty_First_Node   "<< std::endl << Expected_Uncertainty_First_Node_Determinant<< std::endl << std::endl;	

				std::cout << "Last_node_Uncertainty_Determinant   "<< std::endl << Last_node_Uncertainty_Determinant << std::endl;			
				std::cout << "Expected_Uncertainty_Absolute   "<< std::endl << Expected_Uncertainty_Absolute_Determinant << std::endl << std::endl;

				std::cout << "Last_node_probability   "<< Last_node_probability << std::endl << std::endl;		

				std::cout << "D_Optimal_Criteria   "<< D_Optimal_Criteria << std::endl;		
				std::cout << "Work   "<< Work << std::endl;	
				std::cout << "Distance   "<< Distance << std::endl;	
		
		}
	      
};
	



void readNodeProbability(	float value[]){
		

	int i=0;
    std::string STRING;
	std::ifstream infile;
	infile.open ("estimatedProb.txt");
//	infile.open ("EstimatedProbBig.txt");
    while(!infile.eof() ) 
        {
			infile >> value[i];
			value[i]=value[i]/100;
	        i++;
        }
	infile.close();
}


std::vector< std::vector<float> > read_Several_Trayectories(){
//	int i=0;
    std::string name;
    
	std::ifstream infile;
	std::vector< std::vector<float> > List_of_Trajectories;

	float a;
    
    int j=1;
    char char_name[50];

	sprintf (char_name, "paths/path_%d.txt", j);
//	std::cout << "file name:  " << char_name << std::endl;



	while (std::ifstream(char_name))
	{
		std::vector<float> Reader;
		std::cout  <<"Reading "<< char_name << std::endl;
		infile.open (char_name);


		
		while(!infile.eof() ) 
        {
			infile >> a;
			Reader.push_back(a);
//	        std::cerr << a << std::endl;
//			std::cerr << Reader.at(i) << std::endl;
//			i++;
        }
//    std::cerr << Reader << std::endl;
		List_of_Trajectories.push_back(Reader);


		
		infile.close();
		j++;
		sprintf (char_name, "paths/path_%d.txt", j);
//		std::cout  <<"Cycle Over " << std::endl;
			
	}
	j--;

	return List_of_Trajectories;
}






// Linear Transformation
Eigen::Matrix3d CalculateRelativeJacobian1( double PreviousPose[], double CurrentPose[]){
	Eigen::Matrix3d JacobianRelative1;
	float s1=sin(PreviousPose[2]);
	float c1=cos(PreviousPose[2]);
	
	JacobianRelative1(0,0)= -c1;
	JacobianRelative1(0,1)=-s1;
	JacobianRelative1(0,2)=-(CurrentPose[0]- PreviousPose[0])*s1 + (CurrentPose[1]- PreviousPose[1])*c1;
	
	JacobianRelative1(1,0)=s1;
	JacobianRelative1(1,1)=-c1;
	JacobianRelative1(1,2)=-(CurrentPose[0]- PreviousPose[0])*c1 - (CurrentPose[1]- PreviousPose[1])*s1;
	
	JacobianRelative1(2,0)=0;
	JacobianRelative1(2,1)=0;
	JacobianRelative1(2,2)=-1;
	
	
	return JacobianRelative1;
		}

Eigen::Matrix3d CalculateRelativeJacobian2( double PreviousPose[], double CurrentPose[]){
	Eigen::Matrix3d JacobianRelative2;
	float s1=sin(PreviousPose[2]);
	float c1=cos(PreviousPose[2]);
	
	JacobianRelative2(0,0)= c1;
	JacobianRelative2(0,1)=s1;
	JacobianRelative2(0,2)=0;
		
	JacobianRelative2(1,0)=-s1;
	JacobianRelative2(1,1)=c1;
	JacobianRelative2(1,2)=0;
		
	JacobianRelative2(2,0)=0;
	JacobianRelative2(2,1)=0;
	JacobianRelative2(2,2)=1;
	return JacobianRelative2;
	}
	
	

	
Eigen::Matrix3d CalculateJacobianTransform1( double PreviousPose[], double CurrentPose[]){
	Eigen::Matrix3d JacobianTransform1;
	float s1=sin(PreviousPose[2]);
	float c1=cos(PreviousPose[2]);
	JacobianTransform1=Eigen::Matrix3d::Identity();

	JacobianTransform1(0,2)= -(CurrentPose[1]- PreviousPose[1]);	
	JacobianTransform1(1,2)=(CurrentPose[0]- PreviousPose[0]);
	return JacobianTransform1;
	}
	
Eigen::Matrix3d CalculateJacobianTransform2( double PreviousPose[], double CurrentPose[]){
	Eigen::Matrix3d JacobianTransform2;
	float s1=sin(PreviousPose[2]);
	float c1=cos(PreviousPose[2]);

	JacobianTransform2=Eigen::Matrix3d::Identity();
	JacobianTransform2(0,0)= c1;
	JacobianTransform2(0,1)=-s1;
		
	JacobianTransform2(1,0)=s1;
	JacobianTransform2(1,1)=c1;

	return JacobianTransform2;
}	
	



float RotatedMatrix( Eigen::Matrix3d Uncertainty, double AngleDiff){
	Eigen::Matrix3d Transform;
	float s1=sin(AngleDiff);
	float c1=cos(AngleDiff);

	Transform=Eigen::Matrix3d::Identity();
	Transform(0,0)= c1;
	Transform(0,1)=s1;
		
	Transform(1,0)=-s1;
	Transform(1,1)=c1;
	
	Uncertainty=Transform*Uncertainty*Transform.transpose();
	return Uncertainty(1,1);


	return s1;
}





Eigen::Matrix3d calculate_Relative_Uncertainty(int Previous_Node, int Current_Node, g2o::SparseBlockMatrix<Eigen::MatrixXd>& spBlock2 , g2o::SparseOptimizer& optimizer){
	
	// Calculate Relative Uncertainty 	
	Eigen::Matrix3d Covariance1   = *spBlock2.block(Current_Node,Current_Node);
	Eigen::Matrix3d Covariance12  = *spBlock2.block(Current_Node,Previous_Node);
	Eigen::Matrix3d Covariance2   = *spBlock2.block(Previous_Node,Previous_Node);

	//Calculation Jacobian_Relative(Trajectory[NodeNumber], Trajectory[Instant]);
	double PreviousPose[3];
	double CurrentPose[3];

	optimizer.vertex(Current_Node)->getEstimateData(CurrentPose);
	optimizer.vertex(Previous_Node)->getEstimateData(PreviousPose);
	
	
	Eigen::Matrix3d JacobianRelative1 = CalculateRelativeJacobian1(PreviousPose, CurrentPose);
	Eigen::Matrix3d JacobianRelative2 = CalculateRelativeJacobian2(PreviousPose, CurrentPose);

	
	//(J1 J2) |C1    C12| J1T		
	//		  |C12T  C2 | J2T
	
	//*
	Eigen::Matrix3d RelativeUncertainty = (JacobianRelative1*Covariance1)*JacobianRelative1.transpose()
		+ JacobianRelative2*Covariance2*JacobianRelative2.transpose()
		+ JacobianRelative1*Covariance12*JacobianRelative2.transpose()
		+ JacobianRelative2*Covariance12.transpose()*JacobianRelative1.transpose();
	return RelativeUncertainty;
	
}


Eigen::Matrix3d calculate_Transform_Uncertainty(int Previous_Node, int Current_Node, g2o::SparseOptimizer& optimizer, Eigen::Matrix3d Current_Uncertainty, Eigen::Matrix3d CovarianceTransform){
	
	
	
	//Calculation Jacobian_Relative(Trajectory[NodeNumber], Trajectory[Instant]);
	double PreviousPose[3];
	double CurrentPose[3];

	optimizer.vertex(Current_Node)->getEstimateData(CurrentPose);
	optimizer.vertex(Previous_Node)->getEstimateData(PreviousPose);
	
	
	Eigen::Matrix3d JacobianTransform1 = CalculateJacobianTransform1(PreviousPose, CurrentPose);
	Eigen::Matrix3d JacobianTransform2 = CalculateJacobianTransform2(PreviousPose, CurrentPose);

	
	//(J1 J2) |C1    C12| J1T		
	//		  |C12T  C2 | J2T
	
	//*
	Eigen::Matrix3d PropagatedUncertainty = JacobianTransform2*CovarianceTransform*JacobianTransform2.transpose()
					+JacobianTransform1*Current_Uncertainty*JacobianTransform1.transpose();	
						
	return PropagatedUncertainty;
	
}

double distance_between_nodes(int Previous_Node, int Current_Node, g2o::SparseOptimizer& optimizer){
	double PreviousPose[3];
	double CurrentPose[3];
	float value;

	optimizer.vertex(Current_Node)->getEstimateData(CurrentPose);
	optimizer.vertex(Previous_Node)->getEstimateData(PreviousPose);
	
	value = (CurrentPose[0]-PreviousPose[0])*(CurrentPose[0]-PreviousPose[0]) + (CurrentPose[1]-PreviousPose[1])*(CurrentPose[1]-PreviousPose[1]);
	value = sqrt(value);
	return value;
}




Path_Quantization ExpectedValue(	g2o::SparseBlockMatrix<Eigen::MatrixXd>& spBlock2,  int Trajectory[],int Trajectory_Length, g2o::SparseOptimizer& optimizer){

	double PreviousPose[3];
	double CurrentPose[3];

	Eigen::Matrix3d Relative_Uncertainty;
	Eigen::Matrix3d Transform_Uncertainty;
	Eigen::Matrix3d Complete_Uncertainty;
	
	Eigen::Matrix3d Node_Uncertainty;// = spBlock2(Trajectory[0], Trajectory[0]);
//	Eigen::Matrix3d Node_Uncertainty = Eigen::MatrixXd::Zero(3,3);	
	
	float Relocalization_Probability;
	
	//	Initialization of Covariance transform
	Eigen::Matrix3d CovarianceTransform = Eigen::Matrix3d::Identity();	
	CovarianceTransform(0,0)=0.00025;	CovarianceTransform(1,1)=0.00025; 	CovarianceTransform(2,2)=0.000025;

//Read Probability File
	float NodeProbability[8358];
//	readNodeProbability(NodeProbability);
	//*
	for(int i=0 ; i < 8358 ; i ++){
		NodeProbability[i]=0.8;
	}//*/
	
	

//	Vector of Hypothesis	
	std::vector< Hypothesis >  Hypothesis_vector;
	Hypothesis current_hypothesis(2);
	
	
// Relocalization Probability
	Relative_Uncertainty  = calculate_Relative_Uncertainty(Trajectory[0], Trajectory[1], spBlock2 ,  optimizer);	
	Transform_Uncertainty = calculate_Transform_Uncertainty(Trajectory[0], Trajectory[1],  optimizer,  Eigen::MatrixXd::Zero(3,3) , Relative_Uncertainty);
//	Transform_Uncertainty = calculate_Transform_Uncertainty(Trajectory[0], Trajectory[1],  optimizer,  Eigen::MatrixXd::Zero(3,3) , CovarianceTransform);

	Complete_Uncertainty =   Transform_Uncertainty + Relative_Uncertainty; 
	Relocalization_Probability = pow((NodeProbability[1]),0.5)*erf( 10/(sqrt(2)* RotatedMatrix(Complete_Uncertainty, CurrentPose[2]-PreviousPose[2]) ) ) ;  //FIX THIS!!!


// Hypothesis construction
	// Hypothesis "0"
	current_hypothesis.last_loop_closure_    = 0;
	current_hypothesis.probability_          = 1-Relocalization_Probability;
	current_hypothesis.Absolute_Uncertainty_ = calculate_Transform_Uncertainty(Trajectory[0], Trajectory[1], optimizer,  *spBlock2.block(Trajectory[0], Trajectory[0]) , Relative_Uncertainty);
	current_hypothesis.First_Node_Uncertainty_ = calculate_Transform_Uncertainty(Trajectory[0], Trajectory[1], optimizer,  Eigen::MatrixXd::Zero(3,3) , Relative_Uncertainty);
	current_hypothesis.Relative_Uncertainty_ = Relative_Uncertainty;
		
	Hypothesis_vector.push_back(current_hypothesis);
	
	// Hypothesis "1"
	current_hypothesis.last_loop_closure_      =  1;
	current_hypothesis.probability_            = Relocalization_Probability;
	current_hypothesis.Absolute_Uncertainty_   = *spBlock2.block( Trajectory[1], Trajectory[1]) ;
	current_hypothesis.First_Node_Uncertainty_ = calculate_Relative_Uncertainty(Trajectory[0], Trajectory[1], spBlock2 ,  optimizer);
	current_hypothesis.Relative_Uncertainty_   = Eigen::MatrixXd::Zero(3,3);
	
	Hypothesis_vector.push_back(current_hypothesis);
	
	//DEBUG
//	std::cerr << "Hypothesis  0:  "<<  Hypothesis_vector[0].last_loop_closure_	 << "Hypothesis  1:  "<<  Hypothesis_vector[1].last_loop_closure_<<std::endl;
	
	float probability_acumulator;
	
	for (int i=2; i<Trajectory_Length;i++){
		probability_acumulator=0;
		for (int j=0; j<i;j++){
			Relative_Uncertainty  = calculate_Relative_Uncertainty(Trajectory[i-1], Trajectory[i], spBlock2 ,  optimizer);
			
			Hypothesis_vector[j].Absolute_Uncertainty_ =calculate_Transform_Uncertainty(Trajectory[i-1], Trajectory[i], optimizer,  Hypothesis_vector[j].Absolute_Uncertainty_ , Relative_Uncertainty);
			Hypothesis_vector[j].Relative_Uncertainty_ =calculate_Transform_Uncertainty(Trajectory[i-1], Trajectory[i], optimizer,  Hypothesis_vector[j].Relative_Uncertainty_ , Relative_Uncertainty);			
			
			
			// Relocalization Probability
			Transform_Uncertainty  = calculate_Relative_Uncertainty(Trajectory[j], Trajectory[i], spBlock2 ,  optimizer);	

			Complete_Uncertainty =   Hypothesis_vector[j].Relative_Uncertainty_ + Transform_Uncertainty; 
			Relocalization_Probability = pow((NodeProbability[Trajectory[j]]),0.5)*erf( 10/(sqrt(2)* RotatedMatrix(Complete_Uncertainty, CurrentPose[2]-PreviousPose[2]) ) ) ;  //FIX THIS!!!
			
			Relocalization_Probability = sqrt(Relocalization_Probability);
			
			Hypothesis_vector[j].probability_ = Hypothesis_vector[j].probability_ *(1- Relocalization_Probability);
			
			probability_acumulator = probability_acumulator + Hypothesis_vector[j].probability_;
			//Debug
			//std::cout << "Probabilidad de Fallo   "<<1- Relocalization_Probability<< "  at   "<< j << std::endl;
		}

		// Add Hypothesis "i"
		current_hypothesis.last_loop_closure_      = i;
		current_hypothesis.probability_            = 1-probability_acumulator;
		current_hypothesis.Absolute_Uncertainty_   = *spBlock2.block( Trajectory[i], Trajectory[i]) ;
		current_hypothesis.First_Node_Uncertainty_ = calculate_Relative_Uncertainty(Trajectory[0], Trajectory[i], spBlock2 ,  optimizer);
		current_hypothesis.Relative_Uncertainty_   = Eigen::MatrixXd::Zero(3,3);
		
		Hypothesis_vector.push_back(current_hypothesis);
		
		//Debug
//		std::cout << "current_hypothesis.probability_   "<<current_hypothesis.probability_ << std::endl;
	}
	
	
	////  Expected Value Calculation
	Eigen::Matrix3d Expected_Uncertainty_Relative = Eigen::MatrixXd::Zero(3,3);
	Eigen::Matrix3d Expected_Uncertainty_Absolute = Eigen::MatrixXd::Zero(3,3);
	Eigen::Matrix3d Expected_Uncertainty_First_Node = Eigen::MatrixXd::Zero(3,3);

	// Henry
	Eigen::Matrix3d Marginal_Covariance = Eigen::MatrixXd::Zero(3,3);
	float D_Optimal_Criteria = 0;

	//Andrade	
	float Energy_U_k = 0, previous_Energy_U_k = 0;
	float Work=0;
	float Distance=0;


	for (int i=0; i< Hypothesis_vector.size() ; i++){
		Expected_Uncertainty_Relative = Expected_Uncertainty_Relative + Hypothesis_vector[i].probability_*Hypothesis_vector[i].Relative_Uncertainty_;
		Expected_Uncertainty_Absolute = Expected_Uncertainty_Absolute + Hypothesis_vector[i].probability_*Hypothesis_vector[i].Absolute_Uncertainty_;		
		Expected_Uncertainty_First_Node = Expected_Uncertainty_First_Node + Hypothesis_vector[i].probability_*Hypothesis_vector[i].First_Node_Uncertainty_;
		
		Marginal_Covariance = *spBlock2.block( Trajectory[i], Trajectory[i]);
		
		// Henry
		D_Optimal_Criteria = D_Optimal_Criteria + (*spBlock2.block( Trajectory[i], Trajectory[i])  ).determinant() ;
		
		// Andrade
		if (i>1){
			Relative_Uncertainty  = calculate_Relative_Uncertainty(Trajectory[i-1], Trajectory[i], spBlock2 ,  optimizer);
		//		Energy_U_k = 1/( (CovarianceTransform.inverse() + Marginal_Covariance.inverse() ).determinant() );
			Energy_U_k = 1/( (Relative_Uncertainty.inverse() + Marginal_Covariance.inverse() ).determinant() );
			
			if ( Energy_U_k > previous_Energy_U_k){
				Work = Work + ( Energy_U_k - previous_Energy_U_k);
			}
			previous_Energy_U_k = Energy_U_k;
			
			Distance = Distance + distance_between_nodes(Trajectory[i-1], Trajectory[i], optimizer);
			
		}
		
		//Debug
//		std::cout << "Relative_Uncertainty   "<< Relative_Uncertainty.determinant()<<"  Marginal_Covariance  " << Marginal_Covariance.determinant()  << std::endl;			
	/*	std::cout << "Node   "<<   i << std::endl;
		std::cout << "probability_   "<<   Hypothesis_vector[i].probability_ << std::endl;
		std::cout << "Absolute_Uncertainty_  "<< Hypothesis_vector[i].Absolute_Uncertainty_.determinant()<< std::endl;
		std::cout << "Relative_Uncertainty_  "<< Hypothesis_vector[i].Relative_Uncertainty_.determinant()<< std::endl;
		std::cout << "First_Node_Uncertainty_  "<< Hypothesis_vector[i].First_Node_Uncertainty_.determinant()<< std::endl;*/
	}
	

	
	Path_Quantization Results(3);
	
	Results.Last_Node_Uncertainty_to_First = current_hypothesis.First_Node_Uncertainty_;
	Results.Last_Node_Uncertainty_to_First_Determinant = current_hypothesis.First_Node_Uncertainty_.determinant();			

	Results.Expected_Uncertainty_First_Node = Expected_Uncertainty_First_Node;
	Results.Expected_Uncertainty_First_Node_Determinant = Expected_Uncertainty_First_Node.determinant();			


	Results.Expected_Uncertainty_Absolute = Expected_Uncertainty_Absolute;
	Results.Expected_Uncertainty_Absolute_Determinant = Expected_Uncertainty_Absolute.determinant();			

	Results.Last_node_Uncertainty = (	*spBlock2.block( Trajectory[Trajectory_Length-1], Trajectory[Trajectory_Length-1]) );
	Results.Last_node_Uncertainty_Determinant = (	*spBlock2.block( Trajectory[Trajectory_Length-1], Trajectory[Trajectory_Length-1]) ).determinant();			

	Results.Last_node_probability = current_hypothesis.probability_;		
	
	Results.D_Optimal_Criteria = D_Optimal_Criteria;		
	Results.Work = Work;		
	
	Results.Distance = Distance;
	
	
	Results.print_path_properties();
	
	
	/*
	 * O=Odometry_Covariance (n0  n1)
	 * p=Relocalization_Probability(O, n1)
	 * h0 = (1-p, O)
	 * h1 = (p, sigma1)
	 * 
	 * H = h0, h1
	 * 
	 * for i=2:m do
	 *   	prob_acum=0;
	 * 		for every hypothesis
	 * 			Odometry_Covariance(n(i-1), n(i))
	 * 			P = update_covariance(Ph,O, n(i-1), n(i) )
	 * 			p = relocalization_prob(Ph, ni, T?)
	 * 			h = (probh * (1-p), P)
	 * 			prob_acum = prob_acum + prob_h
	 * 		end
	 * 		
	 * 		prob_i = (1-prob_acum)
	 * ///// add new hypothesis
	 * 		hi = (probi, sigma_i)
	 * 
	 * 		H = H unido a hi
	 * 		if i= m
	 * 			P_expected= Suma de p_h*P_h
	 * 			return P_expected
	 * 		end
	 * end
	 * 
	
	*/
	
	


//	return 42;//PreviousPose[0]*100;
	return Results;

	
}






int main(int argc, char** argv)
{

	g2o::EdgeSE2 edge;
	g2o::SparseOptimizer optimizer;
	
	clock_t begin;
	clock_t end;


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
	
	std::vector<Path_Quantization>	Quantization_History;
	
// Reading Several Trajectories	
	std::vector< std::vector<float> > List_of_Trajectories = read_Several_Trayectories();
		

	
//// Reduce size of Covariance Request 
/*
	std::vector <float> covariance_request_list;
	std::vector<float>::iterator it;
	
	covariance_request_list = List_of_Trajectories[0];
	
	for (int i=1;i<List_of_Trajectories.size();i++){
		for(int j=0; j < List_of_Trajectories[i].size();j++){
			it = find (covariance_request_list.begin(), covariance_request_list.end(), List_of_Trajectories[i][j]);
			if (it != covariance_request_list.end())
				covariance_request_list.push_back(List_of_Trajectories[i][j]);
		}
	}
	std::cout<< "covariance_request_list size  " << covariance_request_list.size() << std::endl;
// */	
	for (int k=0;k<List_of_Trajectories.size();k++){
		////Read Trajectory	
		std::vector<float> Readed_Traj;
		
//		readTrajectory(&Readed_Traj);
		Readed_Traj = List_of_Trajectories[k];
		
		std::cout << std::endl<< std::endl<< std::endl<< std::endl<< "Path number "<<k << std::endl;
		
		//Subsampling
		int resample_factor=10;
		std::cout<< "resample_factor  " << resample_factor << std::endl;
		
		div_t divresult;
		divresult = div (Readed_Traj.size(),resample_factor);
		
		int Trayectory_size  = divresult.quot+1;
		int Trajectory[Trayectory_size];
		
		for(int i=0 ; i <  divresult.quot ; i ++){
			Trajectory[i] = Readed_Traj.at(i*resample_factor);
		//		std::cerr << Trajectory[i] << std::endl;
		}
		Trajectory[divresult.quot]= Readed_Traj.back();
		
		
		
		//////////  Retrieval list & Compute Marginals
		
		std::vector< std::pair<int,int> > retrievalList_Big;
		//	int idx_Big[Trayectory_size*Trayectory_size] ;	
		//	int idy_Big[Trayectory_size*Trayectory_size];
		
		std::vector<int> idx_Big;
		idx_Big.resize(Trayectory_size*Trayectory_size);
		
		std::vector<int> idy_Big;
		idy_Big.resize(Trayectory_size*Trayectory_size);
		
		
		for(int i=0 ; i < Trayectory_size ; i ++){
			for (int j=0 ; j < Trayectory_size ; j ++){
				idx_Big[i+Trayectory_size*j] = Trajectory[i];
				idy_Big[i+Trayectory_size*j] = Trajectory[j];
		
			}
		}		
		for(int i=0 ; i < Trayectory_size*Trayectory_size ; i ++)
			retrievalList_Big.push_back( std::pair<int,int>(idx_Big[i],idy_Big[i]));
		
		
		
		begin = clock();
		g2o::SparseBlockMatrix<Eigen::MatrixXd> spinv_Big(idx_Big.data(),idy_Big.data(),Trayectory_size*Trayectory_size , Trayectory_size*Trayectory_size);
		optimizer.solver()->computeMarginals(spinv_Big,retrievalList_Big);
		end = clock();
		std::cout << "Time to compute marginals is "<<double(end - begin) / CLOCKS_PER_SEC << " seconds with  "<<  Trayectory_size << " nodes" << std::endl;
		
		//BIG FUNCTION TO CALCULATE!!!!!!!!!!!!!!!!!!
		
		begin = clock();
		Quantization_History.push_back(    ExpectedValue(spinv_Big, Trajectory, Trayectory_size, optimizer)   );	
//		std::cerr<<   ExpectedValue(spinv_Big, Trajectory, Trayectory_size, optimizer)    << std::endl;
		end = clock();	
		std::cout << "Time to compute the expected uncertainty is "<<double(end - begin) / CLOCKS_PER_SEC << " seconds  with "<<  Trayectory_size << " nodes" << std::endl;
	
	}



}

