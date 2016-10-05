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
	



void readNodeProbability(	float value[]){
		

	int i=0;
    std::string STRING;
	std::ifstream infile;
	infile.open ("estimatedProb.txt");
    while(!infile.eof() ) 
        {
			infile >> value[i];
	        i++;
        }
	infile.close();
}


void readTrajectory (	std::vector<float> *Reader){
		
//     std::cerr << "entre" << std::endl;

	int i=0;
    std::string STRING;
	std::ifstream infile;
	float a;
//	infile.open ("ShortPath.txt");
//	infile.open ("LongPath.txt");
	infile.open ("paths/path_1.txt");
    while(!infile.eof() ) 
        {
			infile >> a;
			Reader->push_back(a);
//	        std::cerr << a << std::endl;
//			std::cerr << Reader.at(i) << std::endl;
			i++;
        }
//    std::cerr << Reader << std::endl;
	infile.close();
	
	//node number for reentry 3375
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
	
	

// Unscented transformed
Eigen::Matrix3d CalculateUncertaintyRelative( double PreviousPose[], double CurrentPose[],Eigen::Matrix3d Covariance1, Eigen::Matrix3d Covariance12, Eigen::Matrix3d Covariance2){
	
	Eigen::Matrix3d CovarianceNavigation = Eigen::MatrixXd::Zero(3,3);;
	
	Eigen::Matrix<double, 6, 6> EntryCovariance;


	Eigen::Matrix<double, 6, 1> EntryPoint;	

// EntryPoint << x1, y1, theta1, x2, y2, theta2

	for(int i=0 ; i < 3 ; i ++){
				EntryPoint(i) = PreviousPose[i];
	}
//	std::cerr<< "Sigma Points   "<< std::endl << SigmaPoints << std::endl;	
	for(int i=3 ; i < 6 ; i ++){
				EntryPoint(i) = CurrentPose[i-3];
	}
//	std::cerr<< "Point to Transform   "<< std::endl << EntryPoint << std::endl;	
	
//	float s1=sin(PreviousPose[2]);
//	float c1=cos(PreviousPose[2]);
		
/*	
	EntryCovariance << 1,2,4,7,11,16,	
					2,13,23,38,58,83,
					4,23,77,122,182,257,
					7,38,122,294,430,600,
					11,58,182,430,855,1180,
					16,83,257,600,1180,2071;
	
	EntryPoint << 2,3,4,5,6,7;
	//*/
	EntryCovariance.block(0,0,3,3)= Covariance1;//Eigen::MatrixXd::Identity(3,3) ;
	EntryCovariance.block(3,0,3,3)= Covariance12.transpose();//Eigen::MatrixXd::Ones(3,3) ;
	EntryCovariance.block(0,3,3,3)= Covariance12;//Eigen::MatrixXd::Random(3,3) ;
	EntryCovariance.block(3,3,3,3)= Covariance2;//5*Eigen::MatrixXd::Identity(3,3) ;
	
	
	
	
	int L=6;//                                 %numer of states


//	std::cerr<< "Entry Covariance   "<< std::endl << EntryCovariance << std::endl;	


	Eigen::Matrix<double, 6, 13> SigmaPoints   = Eigen::MatrixXd::Zero(6,13); //2l+1=13 SIGMA POINTS
	Eigen::Matrix<double, 3, 13> OutputVectors = Eigen::MatrixXd::Zero(3,13);; //2l+1=13 Output Vectors
	
	Eigen::Matrix<double, 3, 1> NewMean = Eigen::MatrixXd::Zero(3,1); //2l+1=13 Output Vectors
	Eigen::Matrix<double, 3, 3> NewCov  = Eigen::MatrixXd::Zero(3,3);; //2l+1=13 Output Vectors

//	SigmaPoints=Eigen::MatrixXd::Zero(6,13);
//	OutputVectors=Eigen::MatrixXd::Zero(3,13);



	float alpha=1e-3; //                                %default, tunable
	float ki=-4;         //                              %default, tunable
	float beta=2;  //                                   %default, tunable
	float lambda=alpha*alpha*(L+ki)-L;//                    %scaling factor
	float c=L+lambda; //                                %scaling factor
	
	float Wm[2*L+1];
	float Wc[2*L+1];//%weights for means %weights for covariance
	
	Wm[0] = lambda/c;
	Wc[0] = lambda/c + (1- alpha*alpha  +   beta);	
	
	for(int i=1 ; i < 2*L+1 ; i ++){
		Wm[i]=0.5/c;
		Wc[i]=0.5/c;
	}
//	std::cerr<< "C   "<< std::endl << c << std::endl;
//	std::cerr<< "Lambda   "<< std::endl << lambda << std::endl;
//	std::cerr<< "Weight Mean   "<< std::endl << Wm[0] << std::endl;
//	std::cerr<< "Weight Cov   " << std::endl << Wm[0] << std::endl;

	c=sqrt(c);
//	float deltaT=1;
	
	Eigen::LLT<Eigen::MatrixXd> lltOfA(EntryCovariance);
	Eigen::MatrixXd Std = lltOfA.matrixL();

//	std::cerr<< "Cholesky Decomposition   "<< std::endl << Std << std::endl;	










	
	SigmaPoints.col(0)=EntryPoint;

//	std::cerr<< "C   "<< std::endl << c << std::endl;		
//	std::cerr<< "Sigma Points   "<< std::endl << SigmaPoints << std::endl;	


	for(int i=1 ; i < L+1 ; i ++){
				SigmaPoints.col(i) = EntryPoint + c*Std.col(i-1);
	}
//	std::cerr<< "Sigma Points   "<< std::endl << SigmaPoints << std::endl;	
	for(int i=L+1 ; i < 2*L+1 ; i ++){
				SigmaPoints.col(i) = EntryPoint - c*Std.col(i-(L+1));
	}
//	std::cerr<< "Sigma Points   "<< std::endl << SigmaPoints << std::endl;	

//Transformation Ecuation

// x1 y1 theta1 x2 y2 theta2

//deltaX= ct1*(x2-x1) + st1(y2-y1)
//deltaY= -st1(x2-x1) + ct1(y2-y1)
//deltatheta=  -theta1 + theta2
	float Wacum=0;
	for(int i=0 ; i < 2*L+1 ; i ++){
	
//		OutputVectors(0,i)=  SigmaPoints(3,i)-SigmaPoints(0,i) ;//ct1*(x2-x1) + st1(y2-y1)
//		OutputVectors(1,i)=  SigmaPoints(4,i)-SigmaPoints(1,i) ;

		
		OutputVectors(0,i)=  cos(SigmaPoints(2,i))*(SigmaPoints(3,i)-SigmaPoints(0,i) )   +  sin(SigmaPoints(2,i))*(SigmaPoints(4,i)-SigmaPoints(1,i) )  ;//ct1*(x2-x1) + st1(y2-y1)
		OutputVectors(1,i)= -sin(SigmaPoints(2,i))*(SigmaPoints(3,i)-SigmaPoints(0,i) )   +  cos(SigmaPoints(2,i))*(SigmaPoints(4,i)-SigmaPoints(1,i) );
		OutputVectors(2,i)=	-SigmaPoints(2,i)+SigmaPoints(5,i);			
		
		NewMean = NewMean + Wm[i]*OutputVectors.col(i);
//		Wacum= Wacum + Wm[i];
//		std::cerr<< "New Mean   "<< i+1 << std::endl << NewMean << std::endl;		
//		std::cerr<< "Output Vectors   "<< i+1 << std::endl << OutputVectors.col(i) << std::endl;		
	}


	Eigen::Matrix<double, 3, 13> MatrixMean = Eigen::MatrixXd::Zero(3,13);; //2l+1=13 Output Vectors
	
	for(int i=0 ; i < 2*L+1 ; i ++){
		MatrixMean.col(i)=OutputVectors.col(i)-NewMean;
	}
//	MatrixMean << NewMean, NewMean, NewMean, NewMean, NewMean, NewMean, NewMean, NewMean, NewMean, NewMean, NewMean, NewMean, NewMean;
	

		
//	std::cerr<< "New Mean   "<< std::endl << NewMean << std::endl;		

//	std::cerr<< "Output Vectors   "<< std::endl << OutputVectors << std::endl;		
//	std::cerr<< "Rows   "<< std::endl << OutputVectors.rows() << std::endl;		
//	std::cerr<< "Columns   "<< std::endl << OutputVectors.cols() << std::endl;	









	Eigen::Matrix<double, 13, 13> DiagonalCovariance = Eigen::MatrixXd::Zero(13,13); //2l+1=13 Output Vectors
	
	for(int i=0 ; i < 2*L+1 ; i ++){
		DiagonalCovariance(i,i)=Wc[i];
	}
	
	CovarianceNavigation = MatrixMean*DiagonalCovariance*MatrixMean.transpose();
	
//	CovarianceNavigation = MatrixMean*Wc.asDiagonal()*MatrixMean.transpose();

//	std::cerr<< "Covariance   "<< std::endl << CovarianceNavigation << std::endl;		
	Std = CovarianceNavigation.llt().matrixL();
//	std::cerr<< "Standard Deviation   "<< std::endl << Std << std::endl;			
	
	/*
	EntryCovariance = Eigen::MatrixXd::Zero(6,6);

	std::cerr<< "CovarianceZero   "<< std::endl << EntryCovariance << std::endl;
	EntryCovariance.block(0,0,3,3)= Covariance1;//Eigen::MatrixXd::Identity(3,3) ;
	EntryCovariance.block(3,0,3,3)= Covariance12.transpose();//Eigen::MatrixXd::Ones(3,3) ;
	EntryCovariance.block(0,3,3,3)= Covariance12;//Eigen::MatrixXd::Random(3,3) ;
	EntryCovariance.block(3,3,3,3)= Covariance2;//5*Eigen::MatrixXd::Identity(3,3) ;
	
	
	
	std::cerr<< "CovariancewithIdent   "<< std::endl << EntryCovariance << std::endl;	
	//*/
	
	return CovarianceNavigation;
	}
	
Eigen::Matrix3d CalculateUncertaintyNavigation( double PreviousPose[], double CurrentPose[],Eigen::Matrix3d Covariance1){
	
	Eigen::Matrix3d CovarianceNavigation = Eigen::MatrixXd::Zero(3,3);	
	Eigen::Matrix<double, 6, 6> EntryCovariance;
	Eigen::Matrix<double, 6, 1> EntryPoint;	

	for(int i=0 ; i < 3 ; i ++){
				EntryPoint(i) = PreviousPose[i];
	}
	// CALCULATE TRANSFORMATION FROM PREVIOUS POSE TO CURRENT POSE (Vx,Vy, W) 
	EntryPoint(3)=  cos(PreviousPose[2])*(CurrentPose[0]-PreviousPose[0] )   +  sin(PreviousPose[2])*(CurrentPose[1]- PreviousPose[1] )  ;//ct1*(x2-x1) + st1(y2-y1)
	EntryPoint(4)= -sin(PreviousPose[2])*(CurrentPose[0]-PreviousPose[0] )   +  cos(PreviousPose[2])*(CurrentPose[1]- PreviousPose[1] );
	EntryPoint(5)=	-PreviousPose[2]+ CurrentPose[2];	


//	for(int i=3 ; i < 6 ; i ++){
//				EntryPoint(i) = CurrentPose[i-3];
//	}
	
	
	Eigen::Matrix3d CovarianceTransform = Eigen::Matrix3d::Identity();	
	CovarianceTransform(0,0)=0.00025;
	CovarianceTransform(1,1)=0.00025;
	CovarianceTransform(2,2)=0.000025;
		
	EntryCovariance.block(0,0,3,3)= Covariance1;//Eigen::MatrixXd::Identity(3,3) ;
	EntryCovariance.block(3,0,3,3)= Eigen::MatrixXd::Zero(3,3) ;
	EntryCovariance.block(0,3,3,3)= Eigen::MatrixXd::Zero(3,3) ;
	EntryCovariance.block(3,3,3,3)= CovarianceTransform;//5*Eigen::MatrixXd::Identity(3,3) ;
	

	
	
	int L=6;//                                 %numer of states

	Eigen::Matrix<double, 6, 13> SigmaPoints   = Eigen::MatrixXd::Zero(6,13); //2l+1=13 SIGMA POINTS
	Eigen::Matrix<double, 3, 13> OutputVectors = Eigen::MatrixXd::Zero(3,13);; //2l+1=13 Output Vectors	
	Eigen::Matrix<double, 3, 1> NewMean = Eigen::MatrixXd::Zero(3,1); //2l+1=13 Output Vectors
	Eigen::Matrix<double, 3, 3> NewCov  = Eigen::MatrixXd::Zero(3,3);; //2l+1=13 Output Vectors

	float alpha=1e-3; //                                %default, tunable
	float ki=-4;         //                              %default, tunable
	float beta=2;  //                                   %default, tunable
	float lambda=alpha*alpha*(L+ki)-L;//                    %scaling factor
	float c=L+lambda; //                                %scaling factor
	
	float Wm[2*L+1];
	float Wc[2*L+1];//%weights for means %weights for covariance
	
	Wm[0] = lambda/c;
	Wc[0] = lambda/c + (1- alpha*alpha  +   beta);	
	
	for(int i=1 ; i < 2*L+1 ; i ++){
		Wm[i]=0.5/c;
		Wc[i]=0.5/c;
	}

	c=sqrt(c);
	Eigen::LLT<Eigen::MatrixXd> lltOfA(EntryCovariance);
	Eigen::MatrixXd Std = lltOfA.matrixL();                                                 //CHOLESKY descomposition

	
	SigmaPoints.col(0)=EntryPoint;
	
	for(int i=1 ; i < L+1 ; i ++){
				SigmaPoints.col(i) = EntryPoint + c*Std.col(i-1);
	}
	for(int i=L+1 ; i < 2*L+1 ; i ++){
				SigmaPoints.col(i) = EntryPoint - c*Std.col(i-(L+1));
	}
//	std::cerr<< "Sigma Points   "<< std::endl << SigmaPoints << std::endl;


//Transformation Ecuation //MODIFY EQUATIO FOR MOVEMENT EQUATIONS!!!!

// x1 y1 theta1 Vx Vy w

//deltaX= ct1*Vx - st1*Vy + x1
//deltaY= st1*Vx + ct1*Vy + y1
//deltatheta=  theta1 + w
//	float Wacum=0;
	for(int i=0 ; i < 2*L+1 ; i ++){
		OutputVectors(0,i)=  cos(SigmaPoints(2,i))*SigmaPoints(3,i)   -  sin(SigmaPoints(2,i))*SigmaPoints(4,i)  + SigmaPoints(0,i) ;//ct1*(x2-x1) + st1(y2-y1)
		OutputVectors(1,i)= sin(SigmaPoints(2,i))*SigmaPoints(3,i)    +  cos(SigmaPoints(2,i))*SigmaPoints(4,i)  + SigmaPoints(1,i);
		OutputVectors(2,i)=	SigmaPoints(2,i)+SigmaPoints(5,i);			
		
		NewMean = NewMean + Wm[i]*OutputVectors.col(i);
	}

//	std::cerr<< "New Mean   "<< std::endl << NewMean << std::endl;		
	Eigen::Matrix<double, 3, 13> MatrixMean = Eigen::MatrixXd::Zero(3,13);; //2l+1=13 Output Vectors
	
	for(int i=0 ; i < 2*L+1 ; i ++){
		MatrixMean.col(i)=OutputVectors.col(i)-NewMean;
	}
	
	Eigen::Matrix<double, 13, 13> DiagonalCovariance = Eigen::MatrixXd::Zero(13,13); //2l+1=13 Output Vectors
	
	for(int i=0 ; i < 2*L+1 ; i ++){
		DiagonalCovariance(i,i)=Wc[i];
	}
	
	CovarianceNavigation = MatrixMean*DiagonalCovariance*MatrixMean.transpose();
//	Std = CovarianceNavigation.llt().matrixL();
	return CovarianceNavigation;
	}

Eigen::Matrix3d CalculateUncertaintyNavigation2( double PreviousPose[], double CurrentPose[],Eigen::Matrix3d CovarianceCurrent,Eigen::Matrix3d CovarianceTransformation){
	
	Eigen::Matrix3d CovarianceNavigation = Eigen::MatrixXd::Zero(3,3);	
	Eigen::Matrix<double, 6, 6> EntryCovariance;
	Eigen::Matrix<double, 6, 1> EntryPoint;	

	for(int i=0 ; i < 3 ; i ++){
				EntryPoint(i) = PreviousPose[i];
	}
	// CALCULATE TRANSFORMATION FROM PREVIOUS POSE TO CURRENT POSE (Vx,Vy, W) 
	EntryPoint(3)=  cos(PreviousPose[2])*(CurrentPose[0]-PreviousPose[0] )   +  sin(PreviousPose[2])*(CurrentPose[1]- PreviousPose[1] )  ;//ct1*(x2-x1) + st1(y2-y1)
	EntryPoint(4)= -sin(PreviousPose[2])*(CurrentPose[0]-PreviousPose[0] )   +  cos(PreviousPose[2])*(CurrentPose[1]- PreviousPose[1] );
	EntryPoint(5)=	-PreviousPose[2]+ CurrentPose[2];	


	EntryCovariance.block(0,0,3,3)= CovarianceCurrent;//Eigen::MatrixXd::Identity(3,3) ;
	EntryCovariance.block(3,0,3,3)= Eigen::MatrixXd::Zero(3,3) ;
	EntryCovariance.block(0,3,3,3)= Eigen::MatrixXd::Zero(3,3) ;
	EntryCovariance.block(3,3,3,3)= CovarianceTransformation;//5*Eigen::MatrixXd::Identity(3,3) ;
	

	
	
	int L=6;//                                 %numer of states

	Eigen::Matrix<double, 6, 13> SigmaPoints   = Eigen::MatrixXd::Zero(6,13); //2l+1=13 SIGMA POINTS
	Eigen::Matrix<double, 3, 13> OutputVectors = Eigen::MatrixXd::Zero(3,13);; //2l+1=13 Output Vectors	
	Eigen::Matrix<double, 3, 1> NewMean = Eigen::MatrixXd::Zero(3,1); //2l+1=13 Output Vectors
	Eigen::Matrix<double, 3, 3> NewCov  = Eigen::MatrixXd::Zero(3,3);; //2l+1=13 Output Vectors

	float alpha=1e-3; //                                %default, tunable
	float ki=-4;         //                              %default, tunable
	float beta=2;  //                                   %default, tunable
	float lambda=alpha*alpha*(L+ki)-L;//                    %scaling factor
	float c=L+lambda; //                                %scaling factor
	
	float Wm[2*L+1];
	float Wc[2*L+1];//%weights for means %weights for covariance
	
	Wm[0] = lambda/c;
	Wc[0] = lambda/c + (1- alpha*alpha  +   beta);	
	
	for(int i=1 ; i < 2*L+1 ; i ++){
		Wm[i]=0.5/c;
		Wc[i]=0.5/c;
	}

	c=sqrt(c);
	Eigen::LLT<Eigen::MatrixXd> lltOfA(EntryCovariance);
	Eigen::MatrixXd Std = lltOfA.matrixL();                                                 //CHOLESKY descomposition

	
	SigmaPoints.col(0)=EntryPoint;
	
	for(int i=1 ; i < L+1 ; i ++){
				SigmaPoints.col(i) = EntryPoint + c*Std.col(i-1);
	}
	for(int i=L+1 ; i < 2*L+1 ; i ++){
				SigmaPoints.col(i) = EntryPoint - c*Std.col(i-(L+1));
	}


	for(int i=0 ; i < 2*L+1 ; i ++){
		OutputVectors(0,i)=  cos(SigmaPoints(2,i))*SigmaPoints(3,i)   -  sin(SigmaPoints(2,i))*SigmaPoints(4,i)  + SigmaPoints(0,i) ;//ct1*(x2-x1) + st1(y2-y1)
		OutputVectors(1,i)= sin(SigmaPoints(2,i))*SigmaPoints(3,i)    +  cos(SigmaPoints(2,i))*SigmaPoints(4,i)  + SigmaPoints(1,i);
		OutputVectors(2,i)=	SigmaPoints(2,i)+SigmaPoints(5,i);			
		
		NewMean = NewMean + Wm[i]*OutputVectors.col(i);
	}

	Eigen::Matrix<double, 3, 13> MatrixMean = Eigen::MatrixXd::Zero(3,13);; //2l+1=13 Output Vectors
	
	for(int i=0 ; i < 2*L+1 ; i ++){
		MatrixMean.col(i)=OutputVectors.col(i)-NewMean;
	}
	
	Eigen::Matrix<double, 13, 13> DiagonalCovariance = Eigen::MatrixXd::Zero(13,13); //2l+1=13 Output Vectors
	
	for(int i=0 ; i < 2*L+1 ; i ++){
		DiagonalCovariance(i,i)=Wc[i];
	}
	
	CovarianceNavigation = MatrixMean*DiagonalCovariance*MatrixMean.transpose();
//	Std = CovarianceNavigation.llt().matrixL();
	return CovarianceNavigation;
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




int calcProb(	g2o::SparseBlockMatrix<Eigen::MatrixXd>& spBlock2, int LastNode, int EntryNode, int ExitNode, int Trajectory[],int Trajectory_Length, 	g2o::SparseOptimizer& optimizer, int TrajectoryLenght){
	
	
	int NodeLocalizationProbability[Trajectory_Length];

	Eigen::Matrix3d RelativeUncertainty;
	Eigen::Matrix3d TransformationUncertainty;
	Eigen::Matrix3d PropagatedUncertainty[Trajectory_Length];
	Eigen::Matrix3d PropagatedInternodeUncertainty[Trajectory_Length];
	Eigen::Matrix3d RelativatedUncertainty[Trajectory_Length];
	Eigen::Matrix3d ExpectedValue[TrajectoryLenght];	

	double NodeInstantaneousProbability[Trajectory_Length];
	double LocalizationProbability[Trajectory_Length];
	double AccumulatedProbability[Trajectory_Length];
	
	double PreviousPose[3];
	double CurrentPose[3];


//Read Probability File
	//	readNodeProbability(NodeProbability);
	float NodeProbability[8358];
	for(int i=0 ; i < 8358 ; i ++){
		NodeProbability[i]=1;
	}


//	Initialization of Covariance transform
	Eigen::Matrix3d CovarianceTransform = Eigen::Matrix3d::Identity();	
	CovarianceTransform(0,0)=0.00025;
	CovarianceTransform(1,1)=0.00025;
	CovarianceTransform(2,2)=0.000025;

	//	Definition for JacobianRelative
	Eigen::Matrix3d JacobianRelative1;
	Eigen::Matrix3d JacobianRelative2;

	Eigen::Matrix3d JacobianTransform1;
	Eigen::Matrix3d JacobianTransform2;
	
	// Definition for Calculus of the Transformated Uncertainty	
	Eigen::Matrix3d Covariance1;//  = *spBlock2.block(Trajectory[r-1],Trajectory[r-1]);
	Eigen::Matrix3d Covariance2;//  = *spBlock2.block(Trajectory[r],Trajectory[r]);
	Eigen::Matrix3d Covariance12;// = *spBlock2.block(Trajectory[r-1],Trajectory[r]);
	Eigen::Matrix3d CompleteUncertainty;	
	
	LocalizationProbability[0]=1;
	PropagatedUncertainty[0]=Eigen::Matrix3d::Zero(); //CHECK
	PropagatedInternodeUncertainty[0]=Eigen::Matrix3d::Zero();
	ExpectedValue[0]=Eigen::Matrix3d::Zero();	

	LocalizationProbability[1]=1;
	PropagatedUncertainty[1]=Eigen::Matrix3d::Zero(); //CHECK
	ExpectedValue[1]=Eigen::Matrix3d::Zero();	
	

	
//*	

	
	double PreviousPoseTemp[3];
	double CurrentPoseTemp[3];
	
	Covariance2  =*spBlock2.block(Trajectory[1],Trajectory[1]);
	Covariance12 =*spBlock2.block(Trajectory[0],Trajectory[1]);
	Covariance1  =*spBlock2.block(Trajectory[0],Trajectory[0]);
	optimizer.vertex(Trajectory[1])->getEstimateData(CurrentPose);
	optimizer.vertex(Trajectory[0])->getEstimateData(PreviousPose);
	
	PropagatedInternodeUncertainty[0]= CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);
	PropagatedUncertainty[0]=Covariance2;
	std::cerr<< PropagatedInternodeUncertainty[0].determinant()<<" After Propagation" <<std::endl;
//	std::cerr<< spBlock2.block(Trajectory[1],Trajectory[1])->determinant()<<" In the node" <<std::endl;


// BIG EXPECTED VALUE CALCULATION

	for (int Instant = 1; Instant< Trajectory_Length; Instant++){
		
		ExpectedValue[Instant]=Eigen::Matrix3d::Zero();
		AccumulatedProbability[Instant]=0;
		
		optimizer.vertex(Trajectory[Instant])->getEstimateData(CurrentPose);
		
		Covariance2  = *spBlock2.block(Trajectory[Instant],Trajectory[Instant]);

		optimizer.vertex(Trajectory[Instant-1])->getEstimateData(PreviousPose);
	

		//Calculation of the Transformation Jacobian
		JacobianTransform1 = CalculateJacobianTransform1(PreviousPose, CurrentPose);
		JacobianTransform2 = CalculateJacobianTransform2(PreviousPose, CurrentPose);	
		
		if (Instant==1){
			//Calculation Jacobian_Relative(Trajectory[NodeNumber], Trajectory[Instant]);
			JacobianRelative1 = CalculateRelativeJacobian1(PreviousPose, CurrentPose);
			JacobianRelative2 = CalculateRelativeJacobian2(PreviousPose, CurrentPose);

			// Calculate Relative Uncertainty 	
			Covariance1   = *spBlock2.block(Trajectory[0],Trajectory[0]);
			Covariance12  = *spBlock2.block(Trajectory[0],Trajectory[Instant]);
			
			Covariance1   =	PropagatedUncertainty[0];
			Covariance12   =	PropagatedUncertainty[0];
	

//DEBUG				
			std::cerr<< PropagatedUncertainty[0].determinant()<<" for first time" <<std::endl;
			std::cerr<< PropagatedInternodeUncertainty[0].determinant()<<" First Internode Uncertainty" <<std::endl;
//			PropagatedUncertainty[0]=PropagatedUncertainty[1];
			}
	
//Loop to calculate every node contribution to the expected value
		for (int NodeNumber=0; NodeNumber< Instant; NodeNumber++){
			//Initialization		
			optimizer.vertex(Trajectory[NodeNumber])->getEstimateData(PreviousPose);
			
			//Calculation Jacobian_Relative(Trajectory[NodeNumber], Trajectory[Instant]);
			JacobianRelative1 = CalculateRelativeJacobian1(PreviousPose, CurrentPose);
			JacobianRelative2 = CalculateRelativeJacobian2(PreviousPose, CurrentPose);

			// Calculate Relative Uncertainty 	
			Covariance1   = *spBlock2.block(Trajectory[NodeNumber],Trajectory[NodeNumber]);
			Covariance12  = *spBlock2.block(Trajectory[NodeNumber],Trajectory[Instant]);
			
	//(J1 J2) |C1    C12| J1T		
	//		  |C12T  C2 | J2T

	//*
			RelativeUncertainty = (JacobianRelative1*Covariance1)*JacobianRelative1.transpose()
				+ JacobianRelative2*Covariance2*JacobianRelative2.transpose()
				+ JacobianRelative1*Covariance12*JacobianRelative2.transpose()
				+ JacobianRelative2*Covariance12.transpose()*JacobianRelative1.transpose();
		//*/
				
//			RelativeUncertainty = CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);
				

			if (NodeNumber==0){
//				PropagatedUncertainty[Instant]= RelativeUncertainty;
				PropagatedUncertainty[Instant]= Covariance2;//HAS TO BE BECAUSE THE PROBABILITY IS ROBOCENTRIC AND THE UNCERTAINTY IS MAP RELATED
				PropagatedInternodeUncertainty[Instant]=Eigen::Matrix3d::Zero();//=CovarianceTransform;
//				std::cerr<<"Relative Uncertainty "<< RelativeUncertainty.determinant() << std::endl;	
//				std::cerr<<"  Propagated Uncertainty "<< PropagatedUncertainty[Instant].determinant() << std::endl;
				}



//			Calculate Next Propagated Uncertainty

			PropagatedUncertainty[NodeNumber]=JacobianTransform2*CovarianceTransform*JacobianTransform2.transpose()
							+JacobianTransform1*PropagatedUncertainty[NodeNumber]*JacobianTransform1.transpose();		

	

			PropagatedInternodeUncertainty[NodeNumber] = CalculateUncertaintyNavigation( PreviousPose, CurrentPose, PropagatedInternodeUncertainty[NodeNumber]);		






			
			CompleteUncertainty =   PropagatedInternodeUncertainty[NodeNumber] + RelativeUncertainty; //FIX THIS!!!


			
			//Calculation with error function!		
//			NodeInstantaneousProbability[NodeNumber] =  1-  1*erf( 10/(sqrt(2)*CompleteUncertainty(0,0)) ) * erf( 10/(sqrt(2)*CompleteUncertainty(1,1)) );// * erf( 3.141592/(2*sqrt(2)*CompleteUncertainty(2,2)) );				

			NodeInstantaneousProbability[NodeNumber] = 1 -  pow((NodeProbability[Instant]*1/100),0.5)*erf( 10/(sqrt(2)* RotatedMatrix(CompleteUncertainty, CurrentPose[2]-PreviousPose[2]) ) ) ;				
//			NodeInstantaneousProbability[NodeNumber] = 1 -  erf( 10/(sqrt(2)* RotatedMatrix(CompleteUncertainty, CurrentPose[2]-PreviousPose[2]) ) ) ;				
			
//			NodeInstantaneousProbability[NodeNumber] = 1 -  sqrt(NodeProbability[Instant]*1/100)*erf( 2/(sqrt(2)*CompleteUncertainty(0,0)) ) * erf( 2/(sqrt(2)*CompleteUncertainty(1,1)) );// * erf( 3.141592/(2*sqrt(2)*CompleteUncertainty(2,2)) );				


			LocalizationProbability[NodeNumber] = LocalizationProbability[NodeNumber]*NodeInstantaneousProbability[NodeNumber];
				
			ExpectedValue[Instant] += PropagatedUncertainty[NodeNumber]*LocalizationProbability[NodeNumber];//CHECK THIS
			AccumulatedProbability[Instant] += LocalizationProbability[NodeNumber];






//DEBUG
			if (NodeNumber==0){
				std::cerr<<"Instant Failure Probability;  "<< NodeInstantaneousProbability[NodeNumber] << ";   Acumulated; " << LocalizationProbability[NodeNumber] << ";   Determinant ;" << CompleteUncertainty.determinant() <<std::endl;	
//				std::cerr<<"Instant Failure Probability;  "<< NodeInstantaneousProbability[NodeNumber] << ";   Acumulated; " << LocalizationProbability[NodeNumber] << ";   Directional Uncertainty ;" << RotatedMatrix(CompleteUncertainty, CurrentPose[2]-PreviousPose[2])  <<std::endl;	
//				std::cerr<<"Graph Relative Uncertainty;  "<< RelativeUncertainty.determinant() <<  ";   Odometry Internode Uncertainty; " << PropagatedInternodeUncertainty[NodeNumber].determinant()   <<std::endl;					
//				std::cerr<<"Odometry Uncertainty; " << PropagatedUncertainty[NodeNumber].determinant() << ";   Graph Absolute Uncertainty; " << Covariance2.determinant()   <<std::endl;					
//				std::cerr<<"Instant Failure Probability;  "<< NodeInstantaneousProbability[NodeNumber] <<std::endl;	
//				std::cerr<<"Propagated Uncertainty;  "<< PropagatedUncertainty[NodeNumber].determinant() <<std::endl;
//				std::cerr<< CompleteUncertainty.determinant() <<std::endl;
//				std::cerr<< RelativeUncertainty.determinant() <<std::endl;
//				std::cerr<< PropagatedInternodeUncertainty[NodeNumber].determinant() - Covariance2.determinant()<<std::endl;
				}
				
			if (NodeNumber== Instant-1){
//				std::cerr<<" Instantaneous  " << NodeInstantaneousProbability[NodeNumber] <<std::endl;	
//				std::cerr<<" Propagated Uncertainty Jacobians  " <<std::endl<< PropagatedUncertainty[NodeNumber] <<std::endl;
//				std::cerr<<" Relative Uncertainty Jacobians  " <<std::endl<< RelativeUncertainty <<std::endl;
//				RelativeUncertainty = CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);
//				PropagatedUncertainty[NodeNumber] = CalculateUncertaintyNavigation( PreviousPose, CurrentPose, Covariance1);
//				std::cerr<<" PreviousPose  " <<std::endl<< PreviousPose[0]<<" "<< PreviousPose[1]<<" "<< PreviousPose[2] <<std::endl;
//				std::cerr<<" CurrentPose  "  <<std::endl<< CurrentPose[0] <<" "<< CurrentPose[1] <<" "<< CurrentPose[2]  <<std::endl;
//				std::cerr<<" Relative Uncertainty SUT  " <<std::endl<< RelativeUncertainty <<std::endl;
//				std::cerr<<" Propagated Uncertainty SUT  " <<std::endl<< RelativeUncertainty <<std::endl;
//				std::cerr<< PropagatedInternodeUncertainty[1].determinant()<<" for Internode and this for Total  " << PropagatedUncertainty[1].determinant() <<std::endl;

				}
		 
		}
//DEBUG 
//	std::cerr<< PropagatedUncertainty[Instant].determinant() << std::endl;	
	
		LocalizationProbability[Instant] = 1-AccumulatedProbability[Instant];
//		ExpectedValue[Instant] += LocalizationProbability[Instant]*PropagatedUncertainty[Instant];
		ExpectedValue[Instant] += LocalizationProbability[Instant]*Covariance2;// HAS TO BE!!!!
//		std::cerr<< Covariance2.determinant() << std::endl;	
	}
	
//DEBUG
	for (int i=0; i<TrajectoryLenght; i++){
//			std::cerr<< PropagatedUncertainty[i].determinant() << std::endl;	
//			std::cerr<<"ExpectedValue   "<< ExpectedValue[i].determinant() << std::endl;	
		}
		
	std::cerr<< "Outside Loop Debug" << std::endl;
	std::cerr<< "AccumulatedProbability  "  << AccumulatedProbability[TrajectoryLenght-1] << std::endl;	
	std::cerr<< "LocalizationProbability  " << LocalizationProbability[TrajectoryLenght-1] << std::endl;

	std::cerr<< "Expected Value   "       <<ExpectedValue[TrajectoryLenght-1].determinant() << std::endl;	
	std::cerr<< "Expected Value D-Opt   " << pow(ExpectedValue[TrajectoryLenght-1].determinant(),0.33333333333333333333333) << std::endl;		

	std::cerr<< "Smallest Uncertainty  "<< PropagatedUncertainty[TrajectoryLenght-1].determinant() << std::endl;	
	std::cerr<< "Smallest Uncertainty  D-Opt  "<< pow(PropagatedUncertainty[TrajectoryLenght-1].determinant(),0.33333333333) << std::endl;	
	std::cerr<< "Bigest Uncertainty   "<<PropagatedUncertainty[0].determinant() << std::endl;	

	std::cerr<< "Bigest Uncertainty D-Opt   "<< pow(PropagatedUncertainty[0].determinant(),0.33333333333333333333333) << std::endl;	
	
	PropagatedUncertainty[0]=PropagatedUncertainty[0]*(LocalizationProbability[0]) + 
												PropagatedUncertainty[TrajectoryLenght-1]*(1-LocalizationProbability[0]);
	std::cerr<< "Pseudo Expected   "<<PropagatedUncertainty[0].determinant() << std::endl;		
	std::cerr<< "Pseudo Expected D-Opt   "<< pow(PropagatedUncertainty[0].determinant(),0.33333333333333333333333) << std::endl;	


	return 42;//PreviousPose[0]*100;

	
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






int calcProb2(	g2o::SparseBlockMatrix<Eigen::MatrixXd>& spBlock2,  int Trajectory[],int Trajectory_Length, g2o::SparseOptimizer& optimizer){

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
	//	readNodeProbability(NodeProbability);
	float NodeProbability[8358];
	for(int i=0 ; i < 8358 ; i ++){
		NodeProbability[i]=1;
	}
	
	

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



	for (int i=0; i< Hypothesis_vector.size() ; i++){
		Expected_Uncertainty_Relative = Expected_Uncertainty_Relative + Hypothesis_vector[i].probability_*Hypothesis_vector[i].Relative_Uncertainty_;
		Expected_Uncertainty_Absolute = Expected_Uncertainty_Absolute + Hypothesis_vector[i].probability_*Hypothesis_vector[i].Absolute_Uncertainty_;		
		Expected_Uncertainty_First_Node = Expected_Uncertainty_First_Node + Hypothesis_vector[i].probability_*Hypothesis_vector[i].First_Node_Uncertainty_;
		
		//Debug
	/*	std::cout << "Node   "<<   i << std::endl;
		std::cout << "probability_   "<<   Hypothesis_vector[i].probability_ << std::endl;
		std::cout << "Absolute_Uncertainty_  "<< Hypothesis_vector[i].Absolute_Uncertainty_.determinant()<< std::endl;
		std::cout << "Relative_Uncertainty_  "<< Hypothesis_vector[i].Relative_Uncertainty_.determinant()<< std::endl;
		std::cout << "First_Node_Uncertainty_  "<< Hypothesis_vector[i].First_Node_Uncertainty_.determinant()<< std::endl;*/
	}
	
	//Debug
//	std::cout << "Trajectory_Length   "<<Trajectory_Length << std::endl;
//	std::cout << "Hypothesis Length   "<< Hypothesis_vector.size() << std::endl;
	
	std::cout << "Last Node Uncertainty First   "<< std::endl << current_hypothesis.First_Node_Uncertainty_.determinant() << std::endl;
	std::cout << "Expected_Uncertainty_First_Node   "<< std::endl << Expected_Uncertainty_First_Node.determinant() << std::endl;	
	std::cout << "Expected_Uncertainty_Absolute   "<< std::endl << Expected_Uncertainty_Absolute.determinant() << std::endl;
	std::cout << "Last node Uncertainty   "<< std::endl << (	*spBlock2.block( Trajectory[Trajectory_Length-1], Trajectory[Trajectory_Length-1]) ).determinant() << std::endl;			
	std::cout << "Last node probability   "<< current_hypothesis.probability_ << std::endl;		
	
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
	
	


	return 42;//PreviousPose[0]*100;

	
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
	
	
////Read Trajectory	
	std::vector<float> Readed_Traj;
	readTrajectory(&Readed_Traj);

	//Subsampling
	int resample_factor=30;
	div_t divresult;
	divresult = div (Readed_Traj.size(),resample_factor);
//	std::cerr << "Resampleo   "<< Readed_Traj.size()  <<"  Quotient  "<< divresult.quot << "   Residuo  "<< divresult.rem << std::endl;		

	int Trayectory_size  = divresult.quot;
	int Trajectory[Trayectory_size];

	for(int i=0 ; i <  divresult.quot ; i ++){
		Trajectory[i] = Readed_Traj.at(i*resample_factor);
//		std::cerr << Trajectory[i] << std::endl;
	}



//////////  Retrieval list & Compute Marginals

	std::vector< std::pair<int,int> > retrievalList_Big;
	int idx_Big[Trayectory_size*Trayectory_size] ;	
	int idy_Big[Trayectory_size*Trayectory_size];

	for(int i=0 ; i < Trayectory_size ; i ++){
		for (int j=0 ; j < Trayectory_size ; j ++){
			idx_Big[i+Trayectory_size*j] = Trajectory[i];
			idy_Big[i+Trayectory_size*j] = Trajectory[j];
		}
	}		
	for(int i=0 ; i < Trayectory_size*Trayectory_size ; i ++)
		retrievalList_Big.push_back( std::pair<int,int>(idx_Big[i],idy_Big[i]));

	g2o::SparseBlockMatrix<Eigen::MatrixXd> spinv_Big(idx_Big,idy_Big,Trayectory_size*Trayectory_size , Trayectory_size*Trayectory_size);
	optimizer.solver()->computeMarginals(spinv_Big,retrievalList_Big);

	
//BIG FUNCTION TO CALCULATE!!!!!!!!!!!!!!!!!!

	clock_t begin = clock();	
//	break_point(optimizer);	
//	std::cerr<<   calcProb(spinv_Big, Trajectory[1], Trajectory[2], Trajectory [Trayectory_size], Trajectory, Trayectory_size, optimizer, Trayectory_size)    << std::endl;
	
	std::cerr<<   calcProb2(spinv_Big, Trajectory, Trayectory_size, optimizer)    << std::endl;
	
	clock_t end = clock();
	
	std::cout << "Tarda en la funcion "<<double(end - begin) / CLOCKS_PER_SEC << "segundos  con "<<  Trayectory_size << " nodos" << std::endl;



}

