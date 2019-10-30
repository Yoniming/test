/**
*	输入：fast.fp文件
*	输出：压缩后的out.stream
*	步骤：①提取fast.fp中的第四行质量分数放入矩阵matrix；
②对matrix进行压缩得到out.stream文件；
**/

#include<iostream>
#include "AC.h"
#include<algorithm>
#include<cstring>

//#define Test_Line_maxsize 160
//#define Test_Total_len 4000000
//#define TT 4

using namespace std;

char **input_hist;      //用来存放解压后的矩阵
//char **input_hist_in;  //真实质量分数数据
char **at_cg;   //存放真实的碱基

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
	//===============================1.0 Encoder definition=========================//
	ac_decoder acd_ny;
	string file_str = argv[1];

	int line_maxsize, max_1,min_1,num_symbol, T;     //定义头文件中的参数
	int total_len = 3000000;
	const string file_out = "./"+file_str + ".fastq.stream";
	const string fasta_out = "./"+file_str + ".fasta";
	const string fastq_out = "./"+file_str + ".fastq";
	//const string file_path = "D:\\基因数据\\fastq文件\\" + file_str + ".fastq";
	ac_decoder_init(&acd_ny, file_out.c_str(), &line_maxsize, &max_1,&min_1,&T);
	num_symbol = max_1 - min_1 + 1;
	input_hist = array_char(total_len , line_maxsize);
	at_cg = array_char(total_len, line_maxsize);

	//===============================1.1 Parameter definition & Initialization=======//
	int model_num = num_symbol*num_symbol*num_symbol * 16;
	ac_model * acm_ny = new ac_model[model_num];
	ac_model mean_model;


	int loop_idx = 0;
	int rows = 0;
	int cols = line_maxsize;
	int adaptive_flag = 1;

	long long int sum_size = 0;
	long long int sum = 0;
	char *quick = new char[line_maxsize + 2];
	char table[256];
	//===============================1.2 Encoder definition=========================//

	int * count_char = new int[total_len];	// 动态分配空间 统计每行实际字符个数
	char *input_actg = new char[line_maxsize + 2];
	char *input_buf = new char[line_maxsize + 2];
	char *num_input_hist = new char[total_len];//存放每行的均值

	memset(table, 0, 256);
	//===============================1.3 Encoder Initialization=====================//
	if (T != 2)
	{
		ac_model_init(&mean_model, num_symbol, NULL, adaptive_flag);
	}
	for (int i = 0; i < model_num; i++)
	{
		ac_model_init(&acm_ny[i], num_symbol, NULL, adaptive_flag);
	}
	//==================================从fasta种提取碱基================================
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;
	//=====================================2 Process==============================//
	FILE* fp_out;
	fp_out = fopen("qual.out", "wb");

	int EOF_flag = 0;
	FILE * fp;
	if (!(fp = fopen(fasta_out.c_str(), "r")))
	{
		cout << "Open Failed!" << endl;
		return -1;
	}
	memset(input_actg, 0, line_maxsize + 2);
	memset(input_buf, 0, line_maxsize + 2);
	int blk_cnt = 0;
	while (!EOF_flag) // loop for each canvas
	{
		char * mid, *tmp;
		for (int i = 0; i < total_len; i++)
		{
			fgets(quick, line_maxsize + 20, fp);
			mid = fgets(input_actg, line_maxsize + 2, fp);
			if (feof(fp))
			{
				EOF_flag = 1;
				break;
			}
			count_char[i] = strlen(input_actg) - 1;
			memcpy(at_cg[i], input_actg, count_char[i]);
			rows++;
		}

//===========================================对ERR1638068_1的量化========================================//
		if (T == 0)
		{
			int pre = 0;
			for (int i = 0; i < rows; i++)
			{
				num_input_hist[i] = ac_decode_symbol(&acd_ny, &mean_model);
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
						//int error = ((Q2 - Q1) > 0) ? (Q2 - Q1) : 0;
						//error_rows[i] += error;
						//int E = min(7, error_rows[i] / 8);
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
					input_hist[i][j] = ac_decode_symbol(&acd_ny, &acm_ny[model_idx]);
				}
			}
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						input_hist[i][j] = input_hist[i][j] + min_1;
					}
				}
				for (int i = 0; i < rows; i++)
				{
					fwrite(input_hist[i], sizeof(char), cols, fp_out);
				}
		}
//
//
//
////=====================================测试read—2所用的context量化====================================
		if (T == 1)
		{
			int pre = 0;
			for (int i = 0; i < rows; i++)
			{
				num_input_hist[i] = ac_decode_symbol(&acd_ny, &mean_model);
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
					input_hist[i][j] = ac_decode_symbol(&acd_ny, &acm_ny[model_idx]);
				}
			}
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						input_hist[i][j] = input_hist[i][j] + min_1;
					}
				}
				for (int i = 0; i < rows; i++)
				{
					fwrite(input_hist[i], sizeof(char), cols, fp_out);
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
			input_hist[i][j] = ac_decode_symbol(&acd_ny, &acm_ny[model_idx]);
		}
			}
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
						input_hist[i][j] = input_hist[i][j] + min_1;
				}
			}
			for (int i = 0; i < rows; i++)
			{
				fwrite(input_hist[i], sizeof(char), cols, fp_out);
			}
		}
////===========================================测SRR870667_1所用的context量化========================================//
//		if (T == 3)
//		{
//			int pre = 0;
//			for (int i = 0; i < rows; i++)
//			{
//				num_input_hist[i] = ac_decode_symbol(&acd_ny, &mean_model[i]);
//			}
//			int i = 0;
//			for (int j = 0; j < cols; j++)
//			{
//				for (int k = 0; k < rows; k++)
//				{
//					if (j % 2 == 0)
//					{
//						i = k;
//					}
//					else
//					{
//						i = rows - k - 1;
//					}
//					int F = num_input_hist[i];
//					int model_idx = 0;
//					if (j == 0)
//					{
//						char cur_base = at_cg[i][j];
//						int J0 = table[cur_base];
//						model_idx = J0*num_symbol + F;
//					}
//					else if (j == 1)
//					{
//						char cur_base = at_cg[i][j];
//						char pre_base_1 = at_cg[i][j - 1];
//						int J0 = table[cur_base];
//						int J1 = table[pre_base_1];
//						int G1 = J0 * 4 + J1;
//						int Q1 = input_hist[i][j - 1];
//						if (Q1 == 0)
//						{
//							model_idx = 0;
//						}
//						else
//						{
//							model_idx = F * 20 + G1 + 5 * num_symbol;
//						}
//					}
//					else if (j == 2)
//					{
//						char cur_base = at_cg[i][j];
//						char pre_base_1 = at_cg[i][j - 1];
//						int J0 = table[cur_base];
//						int J1 = table[pre_base_1];
//						int G1 = J0 * 4 + J1;
//						int Q1 = input_hist[i][j - 1];
//						int Q2 = input_hist[i][j - 2];
//						int A = max(Q1, Q2);
//						A = max(A, F);
//						int D = 0;
//						if (Q1 == Q2)
//						{
//							D = 1;
//						}
//						if (Q1 == 0)
//						{
//							model_idx = 0;
//						}
//						else
//						{
//							model_idx = A*num_symbol * 20 * 2 + F * 20 * 2 + G1 * 2 + D + 5 * num_symbol + num_symbol * 20;
//						}
//					}
//					else
//					{
//						char cur_base = at_cg[i][j];
//						char pre_base_1 = at_cg[i][j - 1];
//						int J0 = table[cur_base];
//						int J1 = table[pre_base_1];
//						int G1 = J0 * 4 + J1;
//						int Q1 = input_hist[i][j - 1];
//						int Q2 = input_hist[i][j - 2];
//						int Q3 = input_hist[i][j - 3];
//						int Q4;
//						if (j == 3)
//						{
//							Q4 = Q3;
//						}
//						else
//						{
//							Q4 = input_hist[i][j - 4];
//						}
//						int A = max(Q1, Q2);
//						A = max(A, F);
//						int B = max(Q3, Q4);
//						B = max(B, F);
//						A = max(A, B);
//						//if (Q2 < Q3)
//						//{
//						//	B = Q3;
//						//}
//						int C = 0;
//						if (Q1 == Q2)
//						{
//							C = 1;
//						}
//						int M = j / 13;
//						int a = cols / 13 + 1;
//						if (Q1 == 0)
//						{
//							model_idx = 0;
//						}
//						else
//						{
//							model_idx = A * num_symbol * 2 * a + F * 2 * a + C*a + M + 5 * num_symbol + num_symbol*num_symbol * 40;
//						}
//					}
//					input_hist[i][j] = ac_decode_symbol(&acd_ny, &acm_ny[model_idx]);
//				}
//			}
//						for (int i = 0; i < rows; i++)
			//{
			//	for (int j = 0; j < cols; j++)
			//	{
			//		input_hist[i][j] = input_hist[i][j] + min_1;
			//	}
			//}
//			for (int i = 0; i < rows; i++)
//			{
//				fwrite(input_hist[i], sizeof(char), cols, fp_out);
//			}
//		}
		//============================ destory=====================================
		//ac_encoder_done(&ace_ny);
		
		rows = 0;
	}
	ac_decoder_done(&acd_ny);

	
	for (int i = 0; i < model_num; i++)
	{
		ac_model_done(&acm_ny[i]);
	}
	//fclose(fp_out);
	//===========================合成FASTQ文件========================================
	char buffer_1[160];
	char buffer_2[160];
	char buffer_3[160];
	rewind(fp);



	FILE* qual_out;
	qual_out = fopen("qual.out", "rb");

	FILE* q_out;
	if (!(q_out = fopen(fastq_out.c_str(), "wb")))
	{
		cout << "Open Failed!" << endl;
		return -1;
	}


	while (fgets(buffer_1, 160, fp) != NULL)
	{
		fgets(buffer_2, 160, fp);
		int length_1 = strlen(buffer_1);
		int length_2 = strlen(buffer_2);
		fgets(buffer_3, length_2, qual_out);
		if (buffer_2[length_2 - 1] == '\n')
		{
			for (int j = 0; j < length_1 - 1; j++)  // length-1：最后一位是回车'\n'
			{
				fprintf(q_out, "%c", buffer_1[j]);
			}
			fprintf(q_out, "\n");
			for (int j = 0; j < length_2 - 1; j++)
			{
				fprintf(q_out, "%c", buffer_2[j]);
			}
			fprintf(q_out, "\n");
			fprintf(q_out, "+");
			fprintf(q_out, "\n");
			for (int j = 0; j < length_2 - 1; j++)
			{
				fprintf(q_out, "%c", buffer_3[j]);
			}
			fprintf(q_out, "\n");
		}
		else
		{
			cout << "End of file. " << endl;
			for (int j = 0; j < length_1; j++)  // length-1：最后一位不是回车'\n'
			{
				fprintf(q_out, "%c", buffer_1[j]);
			}
			fprintf(q_out, "\n");
			for (int j = 0; j < length_2; j++)
			{
				fprintf(q_out, "%c", buffer_2[j]);
			}
			fprintf(q_out, "\n");
			fprintf(q_out, "+");
			fprintf(q_out, "\n");
			for (int j = 0; j < length_2; j++)
			{
				fprintf(q_out, "%c", buffer_3[j]);
			}
		}
	}
	delete[] acm_ny;
	delete[] count_char;
	delete[] num_input_hist;
	delete[] quick;
	freearay(rows, input_hist);
	freearay(rows, at_cg);
	fclose(qual_out);
	fclose(q_out);
	fclose(fp);
	return 0;
}









