#include "kalaman.h"
#include "arm_math.h"
#include "IMUTASK.H"

#include <stdlib.h>
KalmanFilter_t Kalman0;
KalmanFilter_t Kalman1;
KalmanFilter_t Kalman2;
uint16_t sizeof_float, sizeof_double;
void Kalman_Filter_Init(KalmanFilter_t *kf, uint8_t xhatSize, uint8_t uSize, uint8_t zSize)
{
	sizeof_float = sizeof(float);
	sizeof_double = sizeof(double);
	
	
	kf->xhatSize = xhatSize;
	kf->zSize = zSize;
	
	

	memset(kf->xhat_data, 0, sizeof_float * xhatSize);
	Matrix_Init(&kf->xhat, kf->xhatSize, 1, (float *)kf->xhat_data);

	// xhatminus x(k|k-1)

	memset(kf->xhatminus_data, 0, sizeof_float * xhatSize);
	Matrix_Init(&kf->xhatminus, kf->xhatSize, 1, (float *)kf->xhatminus_data);


	memset(kf->z_data, 0, sizeof_float * zSize);
	Matrix_Init(&kf->z, kf->zSize, 1, (float *)kf->z_data);
	

	memset(kf->P_data, 0, sizeof_float * xhatSize * xhatSize);
	Matrix_Init(&kf->P, kf->xhatSize, kf->xhatSize, (float *)kf->P_data);

	// create covariance matrix P(k|k-1)
	
	memset(kf->Pminus_data, 0, sizeof_float * xhatSize * xhatSize);
	Matrix_Init(&kf->Pminus, kf->xhatSize, kf->xhatSize, (float *)kf->Pminus_data);

	// state transition matrix F FT
	
	memset(kf->F_data, 0, sizeof_float * xhatSize * xhatSize);
	memset(kf->FT_data, 0, sizeof_float * xhatSize * xhatSize);
	Matrix_Init(&kf->F, kf->xhatSize, kf->xhatSize, (float *)kf->F_data);
	Matrix_Init(&kf->FT, kf->xhatSize, kf->xhatSize, (float *)kf->FT_data);
	
	
	
	

	memset(kf->H_data, 0, sizeof_float * zSize * xhatSize);
	memset(kf->HT_data, 0, sizeof_float * xhatSize * zSize);
	Matrix_Init(&kf->H, kf->zSize, kf->xhatSize, (float *)kf->H_data);
	Matrix_Init(&kf->HT, kf->xhatSize, kf->zSize, (float *)kf->HT_data);

	// process noise covariance matrix Q

	memset(kf->Q_data, 0, sizeof_float * xhatSize * xhatSize);
	Matrix_Init(&kf->Q, kf->xhatSize, kf->xhatSize, (float *)kf->Q_data);

	// measurement noise covariance matrix R

	memset(kf->R_data, 0, sizeof_float * zSize * zSize);
	Matrix_Init(&kf->R, kf->zSize, kf->zSize, (float *)kf->R_data);

	// kalman gain K

	memset(kf->K_data, 0, sizeof_float * xhatSize * zSize);
	Matrix_Init(&kf->K, kf->xhatSize, kf->zSize, (float *)kf->K_data);


	Matrix_Init(&kf->S, kf->xhatSize, kf->xhatSize, (float *)kf->S_data);
	Matrix_Init(&kf->temp_matrix, kf->xhatSize, kf->xhatSize, (float *)kf->temp_matrix_data);
	Matrix_Init(&kf->temp_matrix1, kf->xhatSize, kf->xhatSize, (float *)kf->temp_matrix_data1);
	Matrix_Init(&kf->temp_vector, kf->xhatSize, 1, (float *)kf->temp_vector_data);
	Matrix_Init(&kf->temp_vector1, kf->xhatSize, 1, (float *)kf->temp_vector_data1);

	
	kf->xhat_data[0]=0.001f;
	kf->xhat_data[1]=0.001f;
	kf->P_data[0]=0.1f;
	kf->P_data[1]=0.1f;
	kf->Q_data[0]=0.01f;
//	kf->Q_data[1]=0.001f;
//	kf->Q_data[2]=0.001f;
//	kf->Q_data[3]=0.001f;
	kf->R_data[0]=1.51f;
	kf->R_data[1]=0;
	kf->R_data[2]=0;
	kf->R_data[3]=1.51f;
	kf->K_data[0]=0.1f;
	kf->K_data[1]=0;
	kf->K_data[2]=0;
	kf->K_data[3]=0.1f;

	
	
}





void Kalaman_feedback(KalmanFilter_t *kf,float dt,float z0,float z1)
{
	
	//卡尔曼滤波器黄金五式实现，多能在陀螺仪等代码中找到，此处使用了二阶参数的卡尔曼滤波器做数据融合，可以看bilibili  dr_can控制类视频了解

	
	kf->F_data[0]=1;
	kf->F_data[1]=dt;
	kf->F_data[2]=0;
	kf->F_data[3]=1;
	kf->H_data[0]=1;
	kf->H_data[1]=0;
	kf->H_data[2]=0;
	kf->H_data[3]=1;
	
	kf->z_data[0]=z0;
	kf->z_data[1]=z1;
				
	
	
	
	
	// 先验估计
	// 1. xhat'(k)= A·xhat(k-1) + B·u
	Matrix_Multiply(&kf->F, &kf->xhat, &kf->xhatminus);
	
	
	
	// 预测更新
  // 2. P'(k) = A·P(k-1)·AT + Q
	Matrix_Transpose(&kf->F, &kf->FT);
	Matrix_Multiply(&kf->F, &kf->P, &kf->Pminus);
	kf->temp_matrix.numRows = kf->Pminus.numRows;
	kf->temp_matrix.numCols = kf->FT.numCols;
	Matrix_Multiply(&kf->Pminus, &kf->FT, &kf->temp_matrix); // temp_matrix = F P(k-1) FT
	Matrix_Add(&kf->temp_matrix, &kf->Q, &kf->Pminus);

	
	
	// 量测更新
  // 3. K(k) = P'(k)·HT / (H·P'(k)·HT + R)
	Matrix_Transpose(&kf->H,&kf->HT);
	kf->temp_matrix.numRows = kf->H.numRows;
  kf->temp_matrix.numCols = kf->Pminus.numCols;
	Matrix_Multiply(&kf->H, &kf->Pminus, &kf->temp_matrix);
	kf->temp_matrix1.numRows = kf->temp_matrix.numRows;
  kf->temp_matrix1.numCols = kf->HT.numCols;
	Matrix_Multiply(&kf->temp_matrix, &kf->HT, &kf->temp_matrix1);
	kf->S.numRows = kf->temp_matrix1.numRows;
  kf->S.numCols = kf->R.numCols;
	Matrix_Add(&kf->temp_matrix1, &kf->R, &kf->S);		// S = H P'(k) HT + R
	Matrix_Inverse(&kf->S, &kf->temp_matrix1);
	Matrix_Multiply(&kf->Pminus, &kf->HT, &kf->temp_matrix);
	Matrix_Multiply(&kf->temp_matrix, &kf->temp_matrix1, &kf->K);

	
	// 融合
	// 4. xhat(k) = xhat'(k) + K(k)·(z(k) - H·xhat'(k))
		
	kf->temp_vector.numRows = kf->H.numRows;
	kf->temp_vector.numCols = 1;
	Matrix_Multiply(&kf->H, &kf->xhatminus, &kf->temp_vector); // temp_vector = H xhat'(k)
	kf->temp_vector1.numRows = kf->z.numRows;
	kf->temp_vector1.numCols = 1;
	Matrix_Subtract(&kf->z, &kf->temp_vector, &kf->temp_vector1); // temp_vector1 = z(k) - H·xhat'(k)
	kf->temp_vector.numRows = kf->K.numRows;
	kf->temp_vector.numCols = 1;
	Matrix_Multiply(&kf->K, &kf->temp_vector1, &kf->temp_vector); // temp_vector = K(k)·(z(k) - H·xhat'(k))
	Matrix_Add(&kf->xhatminus, &kf->temp_vector, &kf->xhat);

	
	
	
	// 修正方差
	// 5. P(k) = (1-K(k)·H)·P'(k) ==> P(k) = P'(k)-K(k)·H·P'(k)	
	kf->temp_matrix.numRows = kf->K.numRows;
	kf->temp_matrix.numCols = kf->H.numCols;
	kf->temp_matrix1.numRows = kf->temp_matrix.numRows;
	kf->temp_matrix1.numCols = kf->Pminus.numCols;
	Matrix_Multiply(&kf->K, &kf->H, &kf->temp_matrix);                 // temp_matrix = K(k)·H
	Matrix_Multiply(&kf->temp_matrix, &kf->Pminus, &kf->temp_matrix1); // temp_matrix1 = K(k)·H·P'(k)
	Matrix_Subtract(&kf->Pminus, &kf->temp_matrix1, &kf->P);

	
	
	


}













