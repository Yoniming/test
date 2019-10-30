/**
*	输入：fast.fp文件
*	输出：压缩后的out.stream
*	步骤：①提取fast.fp中的第四行质量分数放入矩阵matrix；
②对matrix进行压缩得到out.stream文件；
**/

#include<iostream>
#include "AC.h"
#include<algorithm>
#include <cstring>
#include <sstream>

#define Test_Line_maxsize 160
#define Test_Total_len 4000000
//#define TT 4

using namespace std;

//动态创建二维数组
char **array_char(int m, int n)
{
	char **arr = new char *[m];
	for (int i = 0; i < m; ++i)
		arr[i] = new char[n];
	return arr;
}
void freearay(int m, char **arr)
{
	for (int i = 0; i < m; ++i)
	{
		delete[]arr[i];
	}
	delete[] arr;
}
int main(int argc, char* argv[])
{
	int T = atoi(argv[2]);
	//===============================通过测试得到NUM_SYMBOL=======//
	//===============================1.1 Parameter definition & Initialization=======//
	//                D:\\基因数据\\fastq文件\\SRR065390_1.fastq
    //                F:\\基因数据\\ERR174310_1.fastq
	string file_str = argv[1];
	const string file_path = "./"+file_str;
	FILE * test_fp;
	if (!(test_fp = fopen(file_path.c_str(), "r")))
	{
		cout << "Open Failed!" << endl;
		return -1;
	}
	char **input_test;
	int rows_test = 0;
	int cols_test = Test_Line_maxsize;
	char * test;
	//===============================1.2 Encoder definition=========================//
	char test_input_buf[Test_Line_maxsize];
	char *test_quick = new char[Test_Line_maxsize];
	int * test_char = new int[Test_Total_len];
	//===============================1.3 Encoder Initialization=====================//
	input_test = array_char(Test_Total_len, Test_Line_maxsize);
	for (int k = 0; k < Test_Total_len; k++)
	{
		memset(input_test[k], 50, Test_Line_maxsize);
	}
	memset(test_input_buf, 0, Test_Line_maxsize);
	memset(test_char, 0, Test_Total_len);
	//=====================================2 Process==============================//
	for (int i = 0; i < Test_Total_len; i++)
	{
		fgets(test_quick, Test_Line_maxsize, test_fp);
		fgets(test_quick, Test_Line_maxsize, test_fp);
		fgets(test_quick, Test_Line_maxsize, test_fp);
		test = fgets(test_input_buf, Test_Line_maxsize, test_fp);
		if (test == NULL)
		{
			break;
		}
		test_char[i] = strlen(test_input_buf) - 1;
		memcpy(input_test[i], test_input_buf, test_char[i]);
		rows_test++;
	}
	int max_1 = 0, min_1 = 100;
	for (int i = 0; i < rows_test; i++)
	{
		for (int j = 0; j < cols_test; j++)
		{
			if (int(input_test[i][j]) >max_1)
				max_1 = int(input_test[i][j]);
		}
	}
	for (int i = 0; i < rows_test; i++)
	{
		for (int j = 0; j < cols_test; j++)
		{
			if (int(input_test[i][j]) <min_1)
				min_1 = int(input_test[i][j]);
		}
	}
	int num_symbol = max_1 - min_1 + 1;
	int line_maxsize = 0;
	for (int i = 0; i < Test_Total_len; i++)
	{
		if (test_char[i]>line_maxsize)
			line_maxsize = test_char[i];
	}
	fclose(test_fp);
	delete[] test_quick;
	delete[] test_char;
	//===============================1.1 Parameter definition & Initialization=======//
	//#define NUM_SYMBOL num_symbol	
	//#define LINE_MAXSIZE line_maxsize      //每行最大长度
	int total_len = 3000000;
	//#define TOTAL_LEN total_len
	
	char **input_hist;
	char **at_cg;

	int model_num = num_symbol*num_symbol*num_symbol * 16;
	int loop_idx = 0;
	char pFileName[100000];

	int rows = 0;
	int cols = line_maxsize;
	int adaptive_flag = 1;

	long long int sum_size = 0;
	long long int sum = 0;
	char *quick = new char[line_maxsize + 2];
	char table[256];
	//===============================1.2 Encoder definition=========================//

	ac_encoder ace_ny;
	ac_model * mean_model = new ac_model[num_symbol];
	ac_model * acm_ny = new ac_model[model_num];

	int * count_char = new int[total_len];	// 动态分配空间 统计每行实际字符个数
	char * input_buf = new char[line_maxsize + 2];
	char *input_actg = new char[line_maxsize + 2];
	int *dy_mean = new int[total_len];         //存放每行的动态均值
	char *num_input_hist = new char[total_len];//存放每行的均值
	int *error_rows = new int[total_len];    //存放每行的误差
	int *dy_error = new int[total_len];        //存放每行的动态均值


	input_hist = array_char(total_len, line_maxsize);
	at_cg = array_char(total_len, line_maxsize);
	memset(table, 0, 256);
	//===============================1.3 Encoder Initialization=====================//
	for (int i = 0; i < num_symbol; i++)
		ac_model_init(&mean_model[i], num_symbol, NULL, adaptive_flag);

	for (int i = 0; i < model_num; i++)
	{
		ac_model_init(&acm_ny[i], num_symbol, NULL, adaptive_flag);
	}
	FILE * fp;
	if (!(fp = fopen(file_path.c_str(), "r")))
	{
		cout << "Open Failed!" << endl;
		return -1;
	}
	//unsigned char *matrix = new unsigned char[LINE_MAXSIZE*TOTAL_LEN];
	//fread(matrix, LINE_MAXSIZE*TOTAL_LEN, sizeof(char), fp);
	//for (int i = 0; i < LINE_MAXSIZE*TOTAL_LEN; i++)
	//{
	//	cout << matrix[i];
	//}
	memset(input_buf, 0, line_maxsize + 2);
	memset(input_actg, 0, line_maxsize + 2);
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;
	//=====================================2 Process==============================//
	
	stringstream ss;
	ss<<T;
	const string file_out = "/media/lioneh/NYTOSHIBA/CODE/"+ file_str + ss.str() + ".stream";
	ac_encoder_init(&ace_ny, file_out.c_str());
	fwrite(&line_maxsize, sizeof(int), 1, ace_ny.fp);
	fwrite(&max_1, sizeof(int), 1, ace_ny.fp);
	fwrite(&min_1, sizeof(int), 1, ace_ny.fp);
	fwrite(&T, sizeof(int), 1, ace_ny.fp);

	while (!feof(fp)) // loop for each canvas
	{
		//if (loop_idx >= TT)
		//	break;
		//==============================2.1 intialization=============================//
		//sprintf(pFileName, "%d.stream", loop_idx++);
		//ac_encoder_init(&ace_ny, pFileName);
		//fwrite(&line_maxsize, sizeof(int), 1, ace_ny.fp);
		//fread

		/*ac_encoder_init(&ace_ny, file_out.c_str());*/

		//memset(num_input_hist, 0, total_len);
		//memset(error_rows, 0, total_len);
		//memset(count_char, 0, total_len);
		////memset(input_buf, 0, LINE_MAXSIZE);

		//for (int k = 0; k < total_len; k++)
		//{
		//	memset(input_hist[k], -1, line_maxsize);
		//}
		//==============================2.2 IO====================================//
		char * tmp;
		char * mid;

		for (int i = 0; i < total_len; i++)
		{
			fgets(quick, line_maxsize + 20, fp);
			mid = fgets(input_actg, line_maxsize + 2, fp);
			fgets(quick, line_maxsize + 2, fp);
			tmp = fgets(input_buf, line_maxsize + 2, fp);
			//fprintf(test_ret, "input_actg=%d\n", input_actg[0]);	
			if (tmp == NULL)
			{
				break;
			}
			count_char[i] = strlen(input_buf) - 1;
			memcpy(at_cg[i], input_actg, count_char[i]);
			memcpy(input_hist[i], input_buf, count_char[i]);
			rows++;
		}

		//for (int i = 0; i < rows; i++)
		//{
		//	for (int j = 0; j < cols; j++)
		//	{
		//		cout << input_hist[i][j];
		//	}
		//	cout << endl;
		//}
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				if (input_hist[i][j] != -1)
					input_hist[i][j] = input_hist[i][j] - min_1;
			}
		}
		
		//===============================2.3  Process=======================//
		// step1: counting mean
		//=====================计算均值=======================
//		int pre = 0;
//		for (int i = 0; i < rows; i++)
//		{
//			float sum = 0.0;
//			for (int j = 0; j < count_char[i]; j++)
//			{
//				sum += input_hist[i][j];
//			}
//			float tmp = (sum) / count_char[i];
//			num_input_hist[i] = tmp;
///*			if (num_input_hist[i] < (num_symbol - 15))
//				num_input_hist[i] = (num_symbol - 15); */   //用于read
//			if (num_input_hist[i] < (num_symbol - 13))
//				num_input_hist[i] = (num_symbol - 13);
//			else if (num_input_hist[i] < (num_symbol - 11))
//				num_input_hist[i] = (num_symbol - 11);
//			else if (num_input_hist[i] < (num_symbol - 9))
//				num_input_hist[i] = (num_symbol - 9);
//			else if (num_input_hist[i] < (num_symbol - 7))
//				num_input_hist[i] = (num_symbol - 7);
//
//			else if (num_input_hist[i] < (num_symbol - 5))
//				num_input_hist[i] = (num_symbol-5);    //用于ERR1638068_1
//
//
//			//if (num_input_hist[i] < 22)
//			//	num_input_hist[i] = 22;
//			//else if (num_input_hist[i] < 24)
//			//	num_input_hist[i] = 24;
//			//else if (num_input_hist[i] < 26)
//			//	num_input_hist[i] = 26;
//			//else if (num_input_hist[i] < 28)
//			//	num_input_hist[i] = 28;
//			//else if (num_input_hist[i] < 30)
//			//	num_input_hist[i] = 30;
//
//			ac_encode_symbol(&ace_ny, &mean_model[pre], num_input_hist[i]);
//		}

		//{
		//	==============output mean==========================//
		//	FILE * fp;
		//	char name[100];
		//	sprintf(name, "mean_%d.dat", loop_idx - 1);
		//	fp = fopen(name, "wb");
		//	fwrite(num_input_hist, sizeof(char), rows, fp);
		//	fclose(fp);
		//}
		//step2: coding loop
//===========================================对ERR1638068_1的量化========================================//
		if (T == 0)
		{
			int pre = 0;
			for (int i = 0; i < rows; i++)
			{
				float sum = 0.0;
				for (int j = 0; j < count_char[i]; j++)
				{
					sum += input_hist[i][j];
				}
				float tmp = (sum) / count_char[i];
				num_input_hist[i] = tmp;
				if (num_input_hist[i] < (num_symbol - 13))
					num_input_hist[i] = (num_symbol - 13);
				else if (num_input_hist[i] < (num_symbol - 11))
					num_input_hist[i] = (num_symbol - 11);
				else if (num_input_hist[i] < (num_symbol - 9))
					num_input_hist[i] = (num_symbol - 9);
				else if (num_input_hist[i] < (num_symbol - 7))
					num_input_hist[i] = (num_symbol - 7);
				else if (num_input_hist[i] < (num_symbol - 5))
					num_input_hist[i] = (num_symbol - 5);
				ac_encode_symbol(&ace_ny, &mean_model[pre], num_input_hist[i]);
			}
			int i = 0;
			for (int j = 0; j < cols; j++)
			{
				for (int k = 0; k < rows; k++)
				{
					if (j % 2 == 0)
					{
						i = k;
					}
					else
					{
						i = rows - k - 1;
					}
					if (input_hist[i][j] == -1)
					{
						continue;
					}
					int F = num_input_hist[i];
					int model_idx = 0, q_model_idx = 0, g_model_idx = 0;
					if (j == 0)
					{
						char cur_base = at_cg[i][j];
						int J0 = table[cur_base];
						model_idx = J0*num_symbol + F;
					}
					else if (j == 1)
					{
						char cur_base = at_cg[i][j];
						char pre_base_1 = at_cg[i][j - 1];
						int J0 = table[cur_base];
						int J1 = table[pre_base_1];
						int G1 = J0 * 4 + J1;
						int Q1 = input_hist[i][j - 1];
						if (Q1 == 0)
						{
							model_idx = G1*num_symbol + F;
						}
						else
						{
							model_idx = F * 20 + G1 + num_symbol + 20 * num_symbol;
						}
					}
					else
					{
						if (j > 80)
							F = F - 1;
						char cur_base = at_cg[i][j];
						char pre_base_1 = at_cg[i][j - 1];
						char pre_base_2 = at_cg[i][j - 2];
						int J0 = table[cur_base];
						int J1 = table[pre_base_1];
						int J2 = table[pre_base_2];
						int G1 = J0 * 4 + J1;
						int G2 = J0 * 4 * 4 + J1 * 4 + J2;
						int Q1 = input_hist[i][j - 1];
						int Q2 = input_hist[i][j - 2];
						int Q3 = 0, Q4 = 0, Q5 = 0;
						if (j == 2)
						{
							Q3 = Q1;
							Q4 = Q2;
							Q5 = Q1;
						}
						else if (j == 3)
						{
							Q3 = input_hist[i][j - 3];
							Q4 = Q1;
							Q5 = Q2;
						}
						else
						{
							Q3 = input_hist[i][j - 3];
							Q4 = input_hist[i][j - 4];
							Q5 = input_hist[i][j - 5];
						}
						int temp[5];
						temp[0] = Q1;
						temp[1] = Q2;
						temp[2] = Q3;
						temp[3] = Q4;
						temp[4] = Q5;
						for (int m = 0; m < 4; m++)
						{
							for (int n = m; n < 4; n++)
							{
								if (temp[m] < temp[n])
								{
									int x = temp[m];
									temp[m] = temp[n];
									temp[n] = temp[m];
								}
							}
						}
						int A = max(Q1, Q2);
						int B = max(Q3, Q4);
						A = max(A, F);
						B = max(B, F);
						int C = 0;
						if (Q1 == Q2)
						{
							C = 1;
						}
						int D = 0;
						if (Q3 == Q4)
						{
							D = 1;
						}
						int M = j / 13;
						int a = 7;
						//M = min(a, M);
						int error = ((Q2 - Q1) > 0) ? (Q2 - Q1) : 0;
						error_rows[i] += error;
						int E = min(7, error_rows[i] / 8);
						if (Q1 == 0)
						{
							model_idx = F * 84 + G2 + 2 * num_symbol + 40 * num_symbol;
						}
						else if (F < (num_symbol - 5))
						{
							A = temp[0];
							q_model_idx = A *(num_symbol - 5) * 8 + F * 8 + M;
							g_model_idx = G1;
							model_idx = q_model_idx * 20 + g_model_idx + 130 * num_symbol;
						}
						else
						{
							q_model_idx = A * 5 * 5 * 2 * 2 + B * 5 * 2 * 2 + F * 2 * 2 + C * 2 + D;
							g_model_idx = G2;
							model_idx = q_model_idx * 84 + g_model_idx;
						}
					}
						ac_encode_symbol(&ace_ny, &acm_ny[model_idx], input_hist[i][j]);
				}
			}
		}



//=====================================测试read—2所用的context量化====================================
		if (T == 1)
		{

			int pre = 0;
			for (int i = 0; i < rows; i++)
			{
				float sum = 0.0;
				for (int j = 0; j < count_char[i]; j++)
				{
					sum += input_hist[i][j];
				}
				float tmp = (sum) / count_char[i];
				num_input_hist[i] = tmp;
				if (num_input_hist[i] < (num_symbol - 15))
					num_input_hist[i] = (num_symbol - 15);
				else if (num_input_hist[i] < (num_symbol - 13))
					num_input_hist[i] = (num_symbol - 13);
				else if (num_input_hist[i] < (num_symbol - 11))
					num_input_hist[i] = (num_symbol - 11);
				else if (num_input_hist[i] < (num_symbol - 9))
					num_input_hist[i] = (num_symbol - 9);
				else if (num_input_hist[i] < (num_symbol - 7))
					num_input_hist[i] = (num_symbol - 7);
				ac_encode_symbol(&ace_ny, &mean_model[pre], num_input_hist[i]);
			}
			int i = 0;
			for (int j = 0; j < cols; j++)
			{
				for (int k = 0; k < rows; k++)
				{
					if (j % 2 == 0)
					{
						i = k;
					}
					else
					{
						i = rows - k - 1;
					}
					if (input_hist[i][j] == -1)
					{
						continue;
					}
					int F = num_input_hist[i];
					int model_idx = 0, q_model_idx = 0, g_model_idx = 0;
					char cur_base = at_cg[i][j];
					int J0 = table[cur_base];
					if (j == 0)
					{
						model_idx = F * 4 + J0;
					}
					else if (j == 1)
					{

						char pre_base_1 = at_cg[i][j - 1];
						int J1 = table[pre_base_1];
						int G1 = J0 * 4 + J1;
						int Q1 = input_hist[i][j - 1];
						if (Q1 < (num_symbol - 13))
							Q1 = (num_symbol - 13);
						else if (Q1 < (num_symbol - 11))
							Q1 = (num_symbol - 11);
						else if (Q1 < (num_symbol - 9))
							Q1 = (num_symbol - 9);
						else if (Q1 < (num_symbol - 7))
							Q1 = (num_symbol - 7);
						else if (Q1 < (num_symbol - 5))
							Q1 = (num_symbol - 5);
						model_idx = Q1 * num_symbol * 16 + F * 16 + G1;
					}
					else
					{
						if (j > 80)
							F = F - 1;
						char cur_base = at_cg[i][j];
						char pre_base_1 = at_cg[i][j - 1];
						char pre_base_2 = at_cg[i][j - 2];
						int J0 = table[cur_base];
						int J1 = table[pre_base_1];
						int J2 = table[pre_base_2];
						int G1 = J0 * 4 + J1;
						int G2 = J0 * 4 * 4 + J1 * 4 + J2;
						int Q1 = input_hist[i][j - 1];
						int Q2 = input_hist[i][j - 2];
						int Q3 = 0, Q4 = 0, Q5 = 0;
						if (j == 2)
						{
							Q3 = Q1;
							Q4 = Q2;
							Q5 = Q1;
						}
						else if (j == 3)
						{
							Q3 = input_hist[i][j - 3];
							Q4 = Q1;
							Q5 = Q2;
						}
						else
						{
							Q3 = input_hist[i][j - 3];
							Q4 = input_hist[i][j - 4];
							Q5 = input_hist[i][j - 5];
						}
						if (Q3 < (num_symbol - 13))
							Q3 = (num_symbol - 13);
						else if (Q3 < (num_symbol - 11))
							Q3 = (num_symbol - 11);
						else if (Q3 < (num_symbol - 9))
							Q3 = (num_symbol - 9);
						else if (Q3 < (num_symbol - 7))
							Q3 = (num_symbol - 7);
						else if (Q3 < (num_symbol - 5))
							Q3 = (num_symbol - 5);
						int temp[5];
						temp[0] = Q1;
						temp[1] = Q2;
						temp[2] = Q3;
						temp[3] = Q4;
						temp[4] = Q5;
						for (int m = 0; m < 4; m++)
						{
							for (int n = m; n < 4; n++)
							{
								if (temp[m] < temp[n])
								{
									int x = temp[m];
									temp[m] = temp[n];
									temp[n] = temp[m];
								}
							}
						}
						int A = max(Q1, Q2);
						int B = max(Q3, Q4);
						A = max(A, F);
						B = max(B, F);
						int C = 0;
						if (Q1 == Q2)
						{
							C = 1;
						}
						int D = 0;
						if (Q3 == Q4)
						{
							D = 1;
						}
						int M = j / 13;
						int a = cols / 13 + 1;
						if (F < 30)
						{
							A = temp[0];
							q_model_idx = A *num_symbol*a + F*a + M;
							g_model_idx = G1;
							model_idx = q_model_idx * 16 + g_model_idx + num_symbol*num_symbol * 16;
						}
						else
						{
							q_model_idx = A * 5 * 5 * 2 * 2 * a + B * 5 * 2 * 2 * a + F * 2 * 2 * a + C * 2 * a + D*a + M;
							g_model_idx = G1;
							model_idx = q_model_idx * 16 + g_model_idx;
						}
					}
					ac_encode_symbol(&ace_ny, &acm_ny[model_idx], input_hist[i][j]);
				}
			}
		}


//=====================================测Illumina所用的context量化====================================
		if (T == 2)
		{
			int i = 0;
			for (int j = 0; j < cols; j++)
			{
				for (int k = 0; k < rows; k++)
				{
					if (j % 2 == 0)
					{
						i = k;
					}
					else
					{
						i = rows - k - 1;
					}
					int model_idx;
					char cur_base = at_cg[i][j];
					int J0 = table[cur_base];
					if (J0 == 4)
					{
						model_idx = 0;
					}
					else
					{
						if (j == 0)
						{
							model_idx = J0 + 1;
						}
						else if (j == 1)
						{
							char pre_base_1 = at_cg[i][j - 1];
							int J1 = table[pre_base_1];
							int G1 = J0 * 5 + J1;
							int Q1 = input_hist[i][j - 1];
							model_idx = Q1 * 20 + G1 + 5;
						}
						else if (j == 2)
						{
							char pre_base_1 = at_cg[i][j - 1];
							int J1 = table[pre_base_1];
							int G1 = J0 * 5 + J1;
							int Q1 = input_hist[i][j - 1];
							int Q2 = input_hist[i][j - 2];
							model_idx = (Q1*num_symbol + Q2) * 20 + G1 + 5 + num_symbol * 20;
						}
						else
						{
							char pre_base_1 = at_cg[i][j - 1];
							char pre_base_2 = at_cg[i][j - 2];
							int J1 = table[pre_base_1];
							int J2 = table[pre_base_2];
							int G1 = J0 * 5 + J1;
							int G2 = J0 * 4 * 4 + J1 * 4 + J2;
							int Q1 = input_hist[i][j - 1];
							int Q2 = input_hist[i][j - 2];
							int Q3 = input_hist[i][j - 3];
							int Q4 = 0;
							if (j == 3)
							{
								Q4 = Q3;
							}
							else
							{
								int Q4 = input_hist[i][j - 4];
							}
							int A = Q1;
							int B = max(Q2, Q3);
							int C = 0;
							if (Q2 == Q3)
							{
								C = 1;
							}
							if (Q3 == Q4)
							{
								C = 2;
							}
							int M = j / 13;
							int a = cols / 13 + 1;
							model_idx = A * num_symbol * 3 * a * 15 + B * 3 * a * 20 + C*a * 20 + M * 20 + G1 + 5 + num_symbol * 20 + num_symbol*num_symbol * 20;
						}
					}
					ac_encode_symbol(&ace_ny, &acm_ny[model_idx], input_hist[i][j]);
					//if (i == rows - 1)
					//	cout << (int)(input_hist[i][j]) << ":    " << model_idx << "   " << J0 << endl;
				}
			}
		}
//===========================================测SRR870667_1所用的context量化========================================//
		if (T == 3)
		{
			int pre = 0;
			for (int i = 0; i < rows; i++)
			{
				float sum = 0.0;
				for (int j = 0; j < count_char[i]; j++)
				{
					sum += input_hist[i][j];
				}
				float tmp = (sum) / count_char[i];
				num_input_hist[i] = tmp;
				if (num_input_hist[i] < (num_symbol - 13))
					num_input_hist[i] = (num_symbol - 13);
				else if (num_input_hist[i] < (num_symbol - 11))
					num_input_hist[i] = (num_symbol - 11);
				else if (num_input_hist[i] < (num_symbol - 9))
					num_input_hist[i] = (num_symbol - 9);
				else if (num_input_hist[i] < (num_symbol - 7))
					num_input_hist[i] = (num_symbol - 7);
				else if (num_input_hist[i] < (num_symbol - 5))
					num_input_hist[i] = (num_symbol - 5);
				ac_encode_symbol(&ace_ny, &mean_model[pre], num_input_hist[i]);
			}
			int i = 0;
			for (int j = 0; j < cols; j++)
			{
				for (int k = 0; k < rows; k++)
				{
					if (j % 2 == 0)
					{
						i = k;
					}
					else
					{
						i = rows - k - 1;
					}
					int F = num_input_hist[i];
					int model_idx = 0;
					if (j == 0)
					{
						char cur_base = at_cg[i][j];
						int J0 = table[cur_base];
						model_idx = J0*num_symbol + F;
					}
					else if (j == 1)
					{
						char cur_base = at_cg[i][j];
						char pre_base_1 = at_cg[i][j - 1];
						int J0 = table[cur_base];
						int J1 = table[pre_base_1];
						int G1 = J0 * 4 + J1;
						int Q1 = input_hist[i][j - 1];
						if (Q1 == 0)
						{
							model_idx = 0;
						}
						else
						{
							model_idx = F * 20 + G1 + 5 * num_symbol;
						}
					}
					else if (j == 2)
					{
						char cur_base = at_cg[i][j];
						char pre_base_1 = at_cg[i][j - 1];
						int J0 = table[cur_base];
						int J1 = table[pre_base_1];
						int G1 = J0 * 4 + J1;
						int Q1 = input_hist[i][j - 1];
						int Q2 = input_hist[i][j - 2];
						int A = max(Q1, Q2);
						A = max(A, F);
						int D = 0;
						if (Q1 == Q2)
						{
							D = 1;
						}
						if (Q1 == 0)
						{
							model_idx = 0;
						}
						else
						{
							model_idx = A*num_symbol * 20 * 2 + F * 20 * 2 + G1 * 2 + D + 5 * num_symbol + num_symbol * 20;
						}
					}
					else
					{
						char cur_base = at_cg[i][j];
						char pre_base_1 = at_cg[i][j - 1];
						int J0 = table[cur_base];
						int J1 = table[pre_base_1];
						int G1 = J0 * 4 + J1;
						int Q1 = input_hist[i][j - 1];
						int Q2 = input_hist[i][j - 2];
						int Q3 = input_hist[i][j - 3];
						int Q4;
						if (j == 3)
						{
							Q4 = Q3;
						}
						else
						{
							Q4 = input_hist[i][j - 4];
						}
						int A = max(Q1, Q2);
						A = max(A, F);
						int B = max(Q3, Q4);
						B = max(B, F);
						A = max(A, B);
						//if (Q2 < Q3)
						//{
						//	B = Q3;
						//}
						int C = 0;
						if (Q1 == Q2)
						{
							C = 1;
						}
						int M = j / 13;
						int a = cols / 13 + 1;
						if (Q1 == 0)
						{
							model_idx = 0;
						}
						else
						{
							model_idx = A * num_symbol * 2 * a + F * 2 * a + C*a + M + 5 * num_symbol + num_symbol*num_symbol * 40;
						}
					}
					ac_encode_symbol(&ace_ny, &acm_ny[model_idx], input_hist[i][j]);
				}
			}
		}
		//============================ destory=====================================
		//ac_encoder_done(&ace_ny);
		rows = 0;
	}
	ac_encoder_done(&ace_ny);
	FILE* fp2;
	fp2 = fopen(file_out.c_str(), "rb");
	fseek(fp2, 0, SEEK_END);
	sum_size = ftell(fp2);
	cout << sum_size << endl;
	//sum += sum_size;
	for (int i = 0; i < model_num; i++)
	{
		ac_model_done(&acm_ny[i]);
	}
	delete[] acm_ny;
	delete[] count_char;
	delete[] num_input_hist;
	delete[] error_rows;
	delete[] quick;
	freearay(rows, input_hist);
	freearay(rows, at_cg);
	//cout << sum << endl;
	fclose(fp);
	return 0;
}




