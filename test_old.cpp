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






void break_point(g2o::SparseOptimizer& optimizer){
	double d, acum=0;
	double PreviousPose[3];
	double CurrentPose[3];
	
	for (int i = 2627; i< 3375; i++){
		optimizer.vertex(i)->getEstimateData(PreviousPose);
		optimizer.vertex(i+1)->getEstimateData(CurrentPose);
		d = (PreviousPose[0]-CurrentPose[0])*(PreviousPose[0]-CurrentPose[0]) 
			+(PreviousPose[1]-CurrentPose[1])*(PreviousPose[1]-CurrentPose[1]);
		acum+=sqrt(d);
		std::cerr<< "Distance  "<< acum << std::endl;
		
		if (acum > 49.14/2){
			std::cerr<< "Node  "<< i << std::endl;
			acum=0;			
		}
		
	}
}




void readNodeProbability(	float value[2732]){
		

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




int calcProb(	g2o::SparseBlockMatrix<Eigen::MatrixXd>& spBlock2, int LastNode, int EntryNode, int ExitNode, int Trajectory[], 	g2o::SparseOptimizer& optimizer, int TrajectoryLenght){
	
	
	int NodeLocalizationProbability[ExitNode-EntryNode];
	

//		std::cerr<< "So far so good 3 Size: "<<	(sizeof(Trajectory)/sizeof(*Trajectory))<<std::endl;		
	
//	int TrajectoryLenght=ExitNode-EntryNode +1;
//	int TrajectoryLenght=340;
	Eigen::Matrix3d RelativeUncertainty;
	Eigen::Matrix3d TransformationUncertainty;
	Eigen::Matrix3d PropagatedUncertainty[TrajectoryLenght];
	Eigen::Matrix3d PropagatedInternodeUncertainty[TrajectoryLenght];
	Eigen::Matrix3d RelativatedUncertainty[TrajectoryLenght];
	Eigen::Matrix3d ExpectedValue[TrajectoryLenght];	

	double NodeInstantaneousProbability[ExitNode-EntryNode+1];
	double LocalizationProbability[ExitNode-EntryNode+1];
	
		
	double AccumulatedProbability[ExitNode-EntryNode+1];
	

	
	double PreviousPose[3];
	double CurrentPose[3];

//Read Probability File
	float NodeProbability[2732];
	readNodeProbability(NodeProbability);
	

	
	
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
	//Average Distance between nodes
	double d, acum=0, anguloMedio=0;
	
	for (int i = 2; i< TrajectoryLenght; i++){
		optimizer.vertex(Trajectory[i-1])->getEstimateData(PreviousPose);
		optimizer.vertex(Trajectory[i])->getEstimateData(CurrentPose);
		d = (PreviousPose[0]-CurrentPose[0])*(PreviousPose[0]-CurrentPose[0]) 
			+(PreviousPose[1]-CurrentPose[1])*(PreviousPose[1]-CurrentPose[1]);
		acum+=sqrt(d);
//		std::cerr<< "Instant w: "<< std::abs(PreviousPose[2]-CurrentPose[2]) << std::endl;

		anguloMedio = anguloMedio + std::abs(PreviousPose[2]-CurrentPose[2]);
		}
//	std::cerr<< "Average Angle between consecutive nodes: "<< anguloMedio/TrajectoryLenght << std::endl;	
//	std::cerr<< "Average Distance between consecutive nodes: "<< acum/TrajectoryLenght << std::endl;	
	acum=acum/TrajectoryLenght;
	//*/




//Calculate open space uncertainty
//Expected Number of Nodes
	optimizer.vertex(Trajectory[0])->getEstimateData(PreviousPose);
	optimizer.vertex(Trajectory[1])->getEstimateData(CurrentPose);
	d = (PreviousPose[0]-CurrentPose[0])*(PreviousPose[0]-CurrentPose[0]) 
			+(PreviousPose[1]-CurrentPose[1])*(PreviousPose[1]-CurrentPose[1]);
	d=sqrt(d);
	
//DEBUG	
//	std::cerr<< "X Previous: "<< PreviousPose[0]<< "   Y Previous: "<< PreviousPose[1] << std::endl;
//	std::cerr<< "X Current: " << CurrentPose[0] << "   Y Current: " << CurrentPose[1] << std::endl;
//	std::cerr<< "Open Space Distance: "<< d << std::endl;	
	float OpenSpaceNodes = d/acum;
	OpenSpaceNodes =OpenSpaceNodes/10;
	std::cerr<< "Expected Number of Nodes: "<< OpenSpaceNodes << std::endl;	




//Open Space Uncertainty

//	PropagatedUncertainty[0]=*spBlock2.block(Trajectory[27],Trajectory[27]);
//	std::cerr<< PropagatedUncertainty[0].determinant()<<" Other Node  " <<std::endl;

//	PropagatedUncertainty[0]=*spBlock2.block(Trajectory[0],Trajectory[0]);
//	std::cerr<< PropagatedUncertainty[0].determinant()<<" BeforePropagation  " <<std::endl;
	
	double PreviousPoseTemp[3];
	double CurrentPoseTemp[3];
	
	
	for (int i = 0; i< OpenSpaceNodes+1; i++){
		
		PreviousPoseTemp[0]=PreviousPose[0]*(1 - i/OpenSpaceNodes) + CurrentPose[0]*(i/OpenSpaceNodes);
		PreviousPoseTemp[1]=PreviousPose[1]*(1 - i/OpenSpaceNodes) + CurrentPose[1]*(i/OpenSpaceNodes);
		PreviousPoseTemp[2]=PreviousPose[2]*(1 - i/OpenSpaceNodes) + CurrentPose[2]*(i/OpenSpaceNodes);

		CurrentPoseTemp[0]=PreviousPose[0]*(1 - (i+1)/OpenSpaceNodes) + CurrentPose[0]*((i+1)/OpenSpaceNodes);
		CurrentPoseTemp[1]=PreviousPose[1]*(1 - (i+1)/OpenSpaceNodes) + CurrentPose[1]*((i+1)/OpenSpaceNodes);
		CurrentPoseTemp[2]=PreviousPose[2]*(1 - (i+1)/OpenSpaceNodes) + CurrentPose[2]*((i+1)/OpenSpaceNodes);

		
		JacobianTransform1 = CalculateJacobianTransform1(PreviousPoseTemp, CurrentPoseTemp);
		JacobianTransform2 = CalculateJacobianTransform2(PreviousPoseTemp, CurrentPoseTemp);
		
		PropagatedUncertainty[0]=JacobianTransform2* CovarianceTransform      *JacobianTransform2.transpose() 
								+JacobianTransform1* PropagatedUncertainty[0] *JacobianTransform1.transpose();	
																
		PropagatedInternodeUncertainty[0]=JacobianTransform2* CovarianceTransform      *JacobianTransform2.transpose() 
								+JacobianTransform1* PropagatedInternodeUncertainty[0] *JacobianTransform1.transpose();	
//	std::cerr<< PropagatedInternodeUncertainty[0].determinant()<<" Inside Loop" <<std::endl;		
		}
	Covariance2  =*spBlock2.block(3660,3660);
	Covariance12 =*spBlock2.block(2732,3660);
	Covariance1  =*spBlock2.block(2732,2732);
	optimizer.vertex(3660)->getEstimateData(CurrentPose);
	optimizer.vertex(2732)->getEstimateData(PreviousPose);
	
	PropagatedInternodeUncertainty[0]= CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);
	PropagatedUncertainty[0]=Covariance2;
	std::cerr<< PropagatedInternodeUncertainty[0].determinant()<<" After Propagation" <<std::endl;
//	std::cerr<< spBlock2.block(Trajectory[1],Trajectory[1])->determinant()<<" In the node" <<std::endl;


















// BIG EXPECTED VALUE CALCULATION

	for (int Instant = 1; Instant< TrajectoryLenght; Instant++){
		
		ExpectedValue[Instant]=Eigen::Matrix3d::Zero();
		AccumulatedProbability[Instant]=0;
		
		optimizer.vertex(Trajectory[Instant])->getEstimateData(CurrentPose);
		
		Covariance2  = *spBlock2.block(Trajectory[Instant],Trajectory[Instant]);

//		std::cerr<< " Node Covariance  " << Covariance2.determinant() << std::endl;	
		
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
	
/*			RelativeUncertainty = (JacobianRelative1*Covariance1)*JacobianRelative1.transpose()
				+ JacobianRelative2*Covariance2*JacobianRelative2.transpose()
				+ JacobianRelative1*Covariance12*JacobianRelative2.transpose()
				+ JacobianRelative2*Covariance12.transpose()*JacobianRelative1.transpose();*/
	
			
//			RelativeUncertainty = RelativeUncertainty;//THIS IS NOT THE EQUATION!!!!
			
//			PropagatedUncertainty[0]=RelativeUncertainty;
//			PropagatedInternodeUncertainty[0]=RelativeUncertainty;
//			PropagatedInternodeUncertainty[0]=PropagatedUncertainty[0];
			
//DEBUG				
//			PropagatedUncertainty[0]=JacobianTransform2*CovarianceTransform*JacobianTransform2.transpose() *	sqrt( (PreviousPose[0]-CurrentPose[0])*(PreviousPose[0]-CurrentPose[0])) / (acum/TrajectoryLenght); 
			std::cerr<< PropagatedUncertainty[0].determinant()<<" for first time" <<std::endl;
			std::cerr<< PropagatedInternodeUncertainty[0].determinant()<<" First Internode Uncertainty" <<std::endl;
//			PropagatedUncertainty[0]=PropagatedUncertainty[1];
			}
			
			
			
	
	
	//DEBUG
/*	
	if(RelativeUncertainty.trace()<0.000525){
		std::cerr<< "Smaller Trace" << RelativeUncertainty.trace()<<"Node"<< Trajectory[Instant]<< std::endl;	
	}
	if(RelativeUncertainty.determinant()<1.5625e-12){
		std::cerr<< "Smaller Determinant" << RelativeUncertainty.determinant()<<"Node"<< Trajectory[Instant]<< std::endl;	
	}
	std::cerr<< "Percentual Determinant Increase  " << 100*RelativeUncertainty.determinant()/1.5625e-12-100<<"%"<< std::endl;	
	std::cerr<< "Percentual Trace Increase  " << 100*RelativeUncertainty.trace()/0.000525-100<<"%"<< std::endl;		
//*/

//	std::cerr<< "Node" << Instant << std::endl;		
		
		
	
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

//			PropagatedUncertainty[NodeNumber] = CalculateUncertaintyNavigation( PreviousPose, CurrentPose, PropagatedUncertainty[NodeNumber]);				
							

//				std::cerr<<"Instant Failure Probability;  "<< NodeInstantaneousProbability[NodeNumber] << ";   Acumulated; " << LocalizationProbability[NodeNumber] << ";   Determinant ;" << CompleteUncertainty.determinant() <<std::endl;	



//			Calculate Next InternodePropagated Uncertainty

//			PropagatedInternodeUncertainty[NodeNumber]=JacobianTransform2*CovarianceTransform*JacobianTransform2.transpose()
//							+JacobianTransform1*PropagatedInternodeUncertainty[NodeNumber]*JacobianTransform1.transpose();		

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



int calcProb2(	g2o::SparseBlockMatrix<Eigen::MatrixXd>& spBlock2, int LastNode, int EntryNode, int ExitNode, int Trajectory[], 	g2o::SparseOptimizer& optimizer, int TrajectoryLenght){
	
	
	Eigen::Matrix3d RelativeUncertainty = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d JumpUncertainty = Eigen::Matrix3d::Zero();
	
	Eigen::Matrix3d PredictedUncertaintyHypothesis[TrajectoryLenght];// = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d PredictedRelativeUncertaintyHypothesis[TrajectoryLenght];// = Eigen::Matrix3d::Zero();
	
	Eigen::Matrix3d ExpectedValue[TrajectoryLenght];// = Eigen::Matrix3d::Zero();	

	double InstantaneousProbability;
	double HypothesisProbability[TrajectoryLenght];
	HypothesisProbability[0]=1;	
	HypothesisProbability[1]=1;	
	double AccumulatedProbability[TrajectoryLenght];
	AccumulatedProbability[0]=0;
	
	double PreviousPose[3];
	double CurrentPose[3];

//Read Probability File
	float NodeProbability[2732];
	readNodeProbability(NodeProbability);
	Eigen::Matrix3d CovarianceTransform = Eigen::Matrix3d::Identity();	

// Definition for Calculus of the Transformed Uncertainty	
	Eigen::Matrix3d Covariance1;
	Eigen::Matrix3d Covariance2;
	Eigen::Matrix3d Covariance12;

	Eigen::Matrix3d CompleteUncertainty;	

	Covariance2  =*spBlock2.block(3660,3660);//Equivalent to the uncertainty of the jump
	Covariance12 =*spBlock2.block(2732,3660);
	Covariance1  =*spBlock2.block(2732,2732);
	optimizer.vertex(3660)->getEstimateData(CurrentPose);
	optimizer.vertex(2732)->getEstimateData(PreviousPose);	

// Initialize Hypothesis
	PredictedUncertaintyHypothesis[0]=Covariance2;
	PredictedRelativeUncertaintyHypothesis[0] = Eigen::Matrix3d::Zero();
	
//FirstPrediction
	PredictedRelativeUncertaintyHypothesis[0]= CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);
	JumpUncertainty = PredictedRelativeUncertaintyHypothesis[0];
	
//	PredictedUncertaintyHypothesis[0] = CalculateUncertaintyNavigation2( PreviousPose, CurrentPose, 	PredictedUncertaintyHypothesis[0], 	PredictedRelativeUncertaintyHypothesis[0]);
//	CompleteUncertainty = CalculateUncertaintyRelative( CurrentPose, CurrentPose, PredictedRelativeUncertaintyHypothesis[0], Eigen::Matrix3d::Zero(), PredictedRelativeUncertaintyHypothesis[0]);
//	InstantaneousProbability =  1-  pow((NodeProbability[Trajectory[1]]*1/100),0.5)  *  erf( 2/(sqrt(2)*CompleteUncertainty(0,0)) )*erf( 2/(sqrt(2)*CompleteUncertainty(1,1)) );// * erf( 3.141592/(2*sqrt(2)*CompleteUncertainty(2,2)) );				

//	HypothesisProbability[0] = HypothesisProbability[0]*InstantaneousProbability;				
	ExpectedValue[0] = Covariance1;



//	std::cerr<<"Jump Uncertainty   "<< JumpUncertainty.determinant()  << "  Predicted  "<<PredictedUncertaintyHypothesis[0].determinant()<< "  El otro  "<< spBlock2.block(3660,3660)->determinant() <<std::endl;







// BIG EXPECTED VALUE CALCULATION


	for (int Node = 1; Node < TrajectoryLenght; Node++){
		
		ExpectedValue[Node]=Eigen::Matrix3d::Zero();
		AccumulatedProbability[Node]=0;
		
		optimizer.vertex(Trajectory[Node])->getEstimateData(CurrentPose);		
		Covariance2  = *spBlock2.block(Trajectory[Node],Trajectory[Node]);

		PredictedUncertaintyHypothesis[Node]=Covariance2;
		PredictedRelativeUncertaintyHypothesis[Node] = Eigen::Matrix3d::Zero();

//Loop to calculate every node contribution to the expected value
		for (int Hypothesis = 0; Hypothesis< Node; Hypothesis++){
						
			if( Hypothesis == 0 ){
				//FOR HYPOTHESIS ZERO (0) EVERYTHING IS BEFORE AND AFTER THE JUMP
			//Initialization		
				optimizer.vertex(Trajectory[1])->getEstimateData(PreviousPose);//After Jump, Entry Node
			// Calculate Relative Uncertainty 	
				Covariance1   = *spBlock2.block(Trajectory[1],Trajectory[1]);//Entry Node
				Covariance12  = *spBlock2.block(Trajectory[1],Trajectory[Node]);//From Entry Node to Current Node
				
				RelativeUncertainty = CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);//From Entry Node to Current Node
							
			// Calculate Next Propagated Uncertainty
				PredictedRelativeUncertaintyHypothesis[0]= CalculateUncertaintyNavigation2( PreviousPose, CurrentPose, 	JumpUncertainty, 	RelativeUncertainty);//First Jump, Then Uncertainty
				PredictedUncertaintyHypothesis[0] = CalculateUncertaintyNavigation2( PreviousPose, CurrentPose,  *spBlock2.block(3660,3660) , 	RelativeUncertainty);//Jump Propagation, then transformation
			
				if(Node==1){			//Because some problems with initialization
					PredictedUncertaintyHypothesis[0] = *spBlock2.block(3660,3660);
					std::cerr<<"FirstTime   "<< RelativeUncertainty.determinant()  << "  Predicted  "<<PredictedUncertaintyHypothesis[0].determinant()<< "  El otro  "<< spBlock2.block(3660,3660)->determinant() <<std::endl;
					}
			
				optimizer.vertex(Trajectory[Hypothesis])->getEstimateData(PreviousPose);			
			// Calculate Relative Uncertainty 	
				Covariance1   = *spBlock2.block(Trajectory[Hypothesis],Trajectory[Hypothesis]);
				Covariance12  = *spBlock2.block(Trajectory[Hypothesis],Trajectory[Node]);
				
				RelativeUncertainty = CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);
			
			}
			
			
			

			else{
		//Initialization		
				optimizer.vertex(Trajectory[Hypothesis])->getEstimateData(PreviousPose);
		// Calculate Relative Uncertainty 	
				Covariance1   = *spBlock2.block(Trajectory[Hypothesis],Trajectory[Hypothesis]);
				Covariance12  = *spBlock2.block(Trajectory[Hypothesis],Trajectory[Node]);
				
				RelativeUncertainty = CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);
		// Calculate Next Propagated Uncertainty
				PredictedRelativeUncertaintyHypothesis[Hypothesis]= CalculateUncertaintyRelative( PreviousPose, CurrentPose, Covariance1, Covariance12, Covariance2);
				PredictedUncertaintyHypothesis[Hypothesis] = CalculateUncertaintyNavigation2( PreviousPose, CurrentPose, 	Covariance1,
															PredictedRelativeUncertaintyHypothesis[Hypothesis]);
			}
			
			
			
			
			
	// Calculate Transformation Uncertainty			
//			CompleteUncertainty =   PredictedRelativeUncertaintyHypothesis[Hypothesis] + RelativeUncertainty; //FIX THIS!!!
			CompleteUncertainty = CalculateUncertaintyRelative( CurrentPose, CurrentPose, PredictedRelativeUncertaintyHypothesis[Hypothesis], Eigen::Matrix3d::Zero(), RelativeUncertainty);
			
			
			
			
			
		//Calculation with error function!		
//			InstantaneousProbability =  1-  pow((NodeProbability[Trajectory[Node]]*1/100),0.1)  *  erf( 2/(sqrt(2)*CompleteUncertainty(0,0)) )*erf( 2/(sqrt(2)*CompleteUncertainty(1,1)) );// * erf( 3.141592/(2*sqrt(2)*CompleteUncertainty(2,2)) );				
			InstantaneousProbability =  1-    erf( 2/(sqrt(2)*CompleteUncertainty(0,0)) )*erf( 2/(sqrt(2)*CompleteUncertainty(1,1)) );// * erf( 3.141592/(2*sqrt(2)*CompleteUncertainty(2,2)) );				
//			InstantaneousProbability = 1 -  pow((NodeProbability[Trajectory[Node]]*1/100),0.5)  *  erf(10/(sqrt(2)* RotatedMatrix(CompleteUncertainty, CurrentPose[2]-PreviousPose[2]) )) ;				

			HypothesisProbability[Hypothesis] = HypothesisProbability[Hypothesis]*InstantaneousProbability;				
			ExpectedValue[Node] += PredictedUncertaintyHypothesis[Hypothesis]*HypothesisProbability[Hypothesis];
			AccumulatedProbability[Node] += HypothesisProbability[Hypothesis];



//DEBUG
			if (Hypothesis==0){
//				std::cerr<<"Instant Failure Probability;  "<< NodeInstantaneousProbability[NodeNumber] << ";   Acumulated; " << LocalizationProbability[NodeNumber] << ";   Determinant ;" << CompleteUncertainty.determinant() <<std::endl;	
//				std::cerr<<"Instant Failure Probability;  "<< InstantaneousProbability << ";   Acumulated; " << HypothesisProbability[Hypothesis] << ";   Directional Uncertainty ;" << RotatedMatrix(CompleteUncertainty, CurrentPose[2]-PreviousPose[2])  << ";   Determinant ;" << CompleteUncertainty.determinant() <<std::endl;	
				std::cerr<<"Graph Relative Uncertainty;  "<< RelativeUncertainty.determinant() <<  ";   Odometry Internode Uncertainty; " << PredictedUncertaintyHypothesis[Hypothesis].determinant()   <<std::endl;					
//				std::cerr<<"Odometry Uncertainty; " << PredictedUncertaintyHypothesis[Hypothesis].determinant() << ";   Graph Absolute Uncertainty; " << Covariance2.determinant()   <<std::endl;					
//				std::cerr<<"Instant Failure Probability;  "<< NodeInstantaneousProbability[NodeNumber] <<std::endl;	
//				std::cerr<<"Propagated Uncertainty;  "<< PropagatedUncertainty[NodeNumber].determinant() <<std::endl;
//				std::cerr<< CompleteUncertainty.determinant() <<std::endl;
//				std::cerr<< RelativeUncertainty.determinant() <<std::endl;
//				std::cerr<< PropagatedInternodeUncertainty[NodeNumber].determinant() - Covariance2.determinant()<<std::endl;
				}
				
			if (Hypothesis== Node-1){
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
	
		HypothesisProbability[Node] = 1-AccumulatedProbability[Node];
//		ExpectedValue[Instant] += LocalizationProbability[Instant]*PropagatedUncertainty[Instant];
		ExpectedValue[Node] += HypothesisProbability[Node]*Covariance2;// HAS TO BE!!!!
//		std::cerr<< Covariance2.determinant() << std::endl;	
	}
	
//DEBUG
	for (int i=0; i<TrajectoryLenght; i++){
//			std::cerr<<"ExpectedValue   "<< ExpectedValue[i].determinant() << " Predicted " << PredictedUncertaintyHypothesis[i].determinant() << std::endl;
//			std::cerr<< PropagatedUncertainty[i].determinant() << std::endl;	
//			std::cerr<<"ExpectedValue   "<< ExpectedValue[i].determinant() << std::endl;	
		}
		
	std::cerr<< "Outside Loop Debug" << std::endl;
	std::cerr<< "Accumulated Probability  "  << AccumulatedProbability[TrajectoryLenght-1] << std::endl;	
	std::cerr<< "Hypothesis Probability  " << HypothesisProbability[TrajectoryLenght-1] << std::endl;

	std::cerr<< "Expected Value   "       <<ExpectedValue[TrajectoryLenght-1].determinant() << std::endl;	
	std::cerr<< "Expected Value D-Opt   " << pow(ExpectedValue[TrajectoryLenght-1].determinant(),0.33333333333333333333333) << std::endl;		

	std::cerr<< "Smallest Uncertainty  "<< PredictedUncertaintyHypothesis[TrajectoryLenght-1].determinant() << std::endl;	
	std::cerr<< "Smallest Uncertainty  D-Opt  "<< pow(PredictedUncertaintyHypothesis[TrajectoryLenght-1].determinant(),0.33333333333) << std::endl;	
	
	std::cerr<< "Bigest Uncertainty   "<< PredictedUncertaintyHypothesis[0].determinant() << std::endl;	
	std::cerr<< "Bigest Uncertainty D-Opt   "<< pow(PredictedUncertaintyHypothesis[0].determinant(),0.33333333333333333333333) << std::endl;	
	
	PredictedUncertaintyHypothesis[0]=PredictedUncertaintyHypothesis[0]*(HypothesisProbability[0]) + 
												PredictedUncertaintyHypothesis[TrajectoryLenght-1]*(1-HypothesisProbability[0]);
	std::cerr<< "Pseudo Expected   "<<PredictedUncertaintyHypothesis[0].determinant() << std::endl;		
	std::cerr<< "Pseudo Expected D-Opt   "<< pow(PredictedUncertaintyHypothesis[0].determinant(),0.33333333333333333333333) << std::endl;	


	return 42;//PreviousPose[0]*100;

	
}








int main(int argc, char** argv)
{

	g2o::EdgeSE2 edge;
	g2o::SparseOptimizer optimizer;

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
	
	std::vector<float> Readed_Traj;
	readTrajectory(&Readed_Traj);
	

	
//Trajectory 3
//*	
	int Difference = Readed_Traj.size();
	int Trajectory[Difference];

	for(int i=0 ; i <  Readed_Traj.size() ; i ++){
		Trajectory[i] = Readed_Traj.at(i);
		std::cerr << Readed_Traj.at(i) << std::endl;
	}

	
	// */
	
	
	optimizer.optimize(10);


	double PreviousPose[3];
	optimizer.vertex(10)->getEstimateData(PreviousPose);
	std::cerr<< PreviousPose[0] <<std::endl;
	std::cerr<< "Hasta aca bien" <<std::endl;

	int LastNode=2732;
	
	int EntryNode=609;
	int ExitNode=1066;
	

//TRAJECTORY 1 !!!	
/*
	int Difference = ExitNode - EntryNode + 2 ;

	int Trajectory[Difference];

	Trajectory[0]=LastNode;
	
	for(int i=1 ; i < Difference ; i ++){
		Trajectory[i]=EntryNode+i-1;
	}		
//*/

//TRAJECTORY 2 !!!	
/*
	int Difference = LastNode-ExitNode+1;
//	int TrajectoryTemp[Difference];
	int Trajectory[241];

	int k=0;
//	for(int i=0 ; i <  Difference/5 ; i ++){
	for(int i=0 ; i <  40 ; i ++){
		Trajectory[i]=LastNode-5*i;
//		std::cerr<< "Sequence  "<< Trajectory[i] <<std::endl;
		k++;
	}
	
//*
 	for(int i=40 ; i <  241 ; i ++){
		Trajectory[i]= 2200 - 5*i;
//		std::cerr<< "Sequence  "<< Trajectory[i] <<std::endl;
		k++;
	}				

	Difference=k;
//	std::cerr<< "Size  "<< k <<std::endl;

//*/

//DEBUG	

//	std::cerr<< "Diference  "<< Difference <<std::endl;
//	std::cerr<< "First  "<< Trajectory[0] <<std::endl;
//	std::cerr<< "Second  "<< Trajectory[1] <<std::endl;
//	std::cerr<< "Last Member  "<< Trajectory[Difference-1] <<std::endl;
		
	std::vector< std::pair<int,int> > retrievalList_Big;
	
//	std::cerr<< "So far so good"<<std::endl;	
	
	int idx_Big[Difference*Difference] ;	
	int idy_Big[Difference*Difference];
	
//	std::cerr<< "So far so good 2"<<std::endl;	

	
	for(int i=0 ; i < Difference ; i ++){
		for (int j=0 ; j < Difference ; j ++){
			idx_Big[i+Difference*j] = Trajectory[i];
			idy_Big[i+Difference*j] = Trajectory[j];

//DEBUG			
//			std::cerr<<idx_Big[i+Difference*j]<<" ";
//			std::cerr<<idy_Big[i+Difference*j]<<std::endl;
		}
	}		
	

		
	for(int i=0 ; i < Difference*Difference ; i ++)
		retrievalList_Big.push_back( std::pair<int,int>(idx_Big[i],idy_Big[i]));

	retrievalList_Big.push_back( std::pair<int,int>(3660,3660) ); //Better than propagate error...make
	retrievalList_Big.push_back( std::pair<int,int>(2732,3660) );

		
//DEBUG		
//	std::cerr<< "So far so good"<<std::endl;


	g2o::SparseBlockMatrix<Eigen::MatrixXd> spinv_Big(idx_Big,idy_Big,Difference*Difference,Difference*Difference);

	std::cerr<< "Hasta aca bien2" <<std::endl;
	optimizer.solver()->computeMarginals(spinv_Big,retrievalList_Big);

//* 
  //DEBUG
	for(int i=0; i < 4 ; i++)
	{
//		std::cerr<<spinv_Big.block(Trajectory[i],Trajectory[i])->determinant()<<std::endl;
//		std::cerr<<spinv.block(i,i)->determinant()<<std::endl;
	}
/*	std::cerr<<spinv.block(1,2)->trace()<<std::endl;

	std::cerr<<spinv_Big.block(2700,645)->trace()<<std::endl;
	std::cerr<<spinv_Big.block(645,2700)->trace()<<std::endl;
*/

//	Eigen::MatrixXd pepito=Eigen::Matrix3d::Zero();




	double pose[3];
	Eigen::Matrix3d MatrixTemplate = Eigen::MatrixXd::Random(3,3);

	
//	optimizer.vertex(10)->getEstimateData(pose);
	
//	CalculateUncertaintyRelative( pose, pose, MatrixTemplate, MatrixTemplate, MatrixTemplate);
	
	
//BIG FUNCTION TO CALCULATE!!!!!!!!!!!!!!!!!!


		std::cerr<< "Hasta aca bien3" <<std::endl;
		
		clock_t begin = clock();
	
	break_point(optimizer);	
//	std::cerr<<   calcProb(spinv_Big, LastNode, EntryNode, ExitNode, Trajectory, optimizer, Difference)    << std::endl;
	
	clock_t end = clock();
	
	std::cout << "Tarda en la funcion "<<double(end - begin) / CLOCKS_PER_SEC << "segundos  con "<<  Difference << " nodos" << std::endl;
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
//	std::cerr<<   pose[0]	   << std::endl;
//	std::cerr<< spinv_Big.block(2700,645)->determinant()<<std::endl;
//	std::cerr<< pepito<<std::endl;
//	std::cerr<< sin(3.141592/4)<<std::endl;
//	std::cerr<< sqrt(3)<<std::endl;	




/*

Eigen::MatrixXd A(3,3);
A << 4,-1,2, -1,6,0, 2,0,5;
std::cerr << "The matrix A is" << std::endl << A << std::endl;
Eigen::LLT<Eigen::MatrixXd> lltOfA(A); // compute the Cholesky decomposition of A
Eigen::MatrixXd L = lltOfA.matrixL(); // retrieve factor L in the decomposition
// The previous two lines can also be written as "L = A.llt().matrixL()"
std::cerr << "The Cholesky factor L is" << std::endl << L << std::endl;
std::cerr << "To check this, let us compute L * L.transpose()" << std::endl;
std::cerr << L * L.transpose() << std::endl;
std::cerr << "This should equal the matrix A" << std::endl;

L = A.llt().matrixL();
std::cerr << "The Cholesky factor L calculated in the other way is" << std::endl << L << std::endl;

//*/




	/* Jacobian from X1 to X2
			x1 		y1 			theta1 					x2 		y2	 theta2
dFx			-cos	-sin	-(x2-x1)sin+(y2-y1)cos		cos		sin		0
dFy			sin		-cos	-(x2-x1)cos+(y2-y1)sin		-sin	cos		0
dFtheta  	0		0			-1						0		0		1

Jacobian from X1 with deltaV and delta Theta
TO CHECK
			x1 		y1 			theta1 					ux 		uy	 utheta
dVx			1		0			-(y2-y1)				cos		-sin	0
dVy			0		1			 x2-x1					sin		cos		0
dVtheta  	0		0				1					0		0		1

*/

	
	
	
	
//	optimizer.save("out.g2o");


}

