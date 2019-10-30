/**
*	输入：fast.fp文件(symbol数量超过40，行数长度超过160且不对齐)
*	输出：压缩后的out.stream
*	步骤：①提取fast.fp中的第四行质量分数放入矩阵matrix；
②对matrix进行压缩得到out.stream文件；
**/


#include<iostream>
#include "AC.h"
#include<algorithm>
#include<vector>
#include<cstring>

using namespace std;

int total_len = 100000;
char **input_hist;
char **at_cg;
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

int de_code(string file_str)
{
	ac_decoder acd_ny;
	int line_max, max_0, min_0, num_symbol, allLen;
	int total_len = 100000;
	//const string file_path = "D:\\基因数据\\fastq文件\\" + file_str + ".fastq";
	const string file_out = file_str + ".fastq.stream";
	const string fasta_out = file_str + ".fasta";
	const string fastq_out = file_str + ".fastq";
	ac_decoder_init(&acd_ny, file_out.c_str(), &line_max, &max_0, &min_0, &allLen);
	int * test_char = new int[allLen-1];
	num_symbol = max_0 - min_0 + 1;
	input_hist = array_char(total_len, line_max);
	at_cg = array_char(total_len, line_max);
	int model_num = num_symbol*num_symbol*num_symbol * 16;
	char *quick = new char[line_max+2];
	char *input_actg = new char[line_max + 2];
	char *input_buf = new char[line_max + 2];
	ac_model * acm_ny = new ac_model[model_num];
	ac_model line_model;
	int loop = 0;
	int rows = 0;
	int adaptive_flag = 1;
	char table[256];
	memset(table, 0, 256);
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;

	ac_model_init(&line_model, line_max+1, NULL, adaptive_flag);
	for (int i = 0; i < model_num; i++)
	{
		ac_model_init(&acm_ny[i], num_symbol, NULL, adaptive_flag);
	}
	FILE* fp_out;
	fp_out = fopen("qual.out", "wb");
	int EOF_flag = 0;
	FILE * fp;
	if (!(fp = fopen(fasta_out.c_str(), "r")))
	{
		cout << "Open Failed!" << endl;
		return -1;
	}
	memset(input_actg, 0, line_max + 2);
	memset(input_buf, 0, line_max + 2);
	for (int i = 0; i < allLen-1; i++)
	{
		test_char[i] = ac_decode_symbol(&acd_ny, &line_model);
	}
	while (!EOF_flag)
	{
		char * mid;
		for (int i = 0; i < total_len; i++)
		{
			fgets(quick, line_max, fp);
			mid = fgets(input_actg, line_max + 2, fp);
			if (feof(fp))
			{
				EOF_flag = 1;
				break;
			}
			memcpy(at_cg[i], input_actg, line_max);
			rows++;
		}

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < test_char[loop*total_len + i]; j++)
				{
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
							//int M = j / 13;
							//int a = test_char[loop*total_len + i] / 13 + 1;
							model_idx = A * num_symbol * 3 * 15 + B * 3 * 20 + C * 20 + G1 + 5 + num_symbol * 20 + num_symbol*num_symbol * 20;
						}
					}
					input_hist[i][j] = ac_decode_symbol(&acd_ny, &acm_ny[model_idx]);
				}
			}
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < test_char[loop*total_len + i]; j++)
				{
					input_hist[i][j] = input_hist[i][j] + min_0;
				}
			}
			for (int i = 0; i < rows; i++)
			{
				fwrite(input_hist[i], sizeof(char), test_char[loop*total_len + i], fp_out);
			}
		rows = 0;
		loop++;
	}
	ac_decoder_done(&acd_ny);
	for (int i = 0; i < model_num; i++)
	{
		ac_model_done(&acm_ny[i]);
	}



	//char buffer_1[2000];
	//char buffer_2[2000];
	//char buffer_3[2000];
	//FILE* qual_out;
	//qual_out = fopen("qual.out", "rb");

	//FILE* q_out;
	//if (!(q_out = fopen(fastq_out.c_str(), "wb")))
	//{
	//	cout << "Open Failed!" << endl;
	//	return -1;
	//}
	//int i = 0;
	//while (fgets(buffer_1, 2000, fp) != NULL)
	//{
	//	fgets(buffer_2, 160, fp);
	//	int length_1 = strlen(buffer_1);
	//	fgets(buffer_2, test_char[i], qual_out);
	//	if (buffer_1[length_1 - 1] == '\n')
	//	{
	//		for (int j = 0; j < length_1 - 1; j++)  // length-1：最后一位是回车'\n'
	//		{
	//			fprintf(q_out, "%c", buffer_1[j]);
	//		}
	//		fprintf(q_out, "\n");
	//		for (int j = 0; j < test_char[i]; j++)
	//		{
	//			fprintf(q_out, "%c", buffer_2[j]);
	//		}
	//		fprintf(q_out, "\n");
	//		fprintf(q_out, "+");
	//		fprintf(q_out, "\n");
	//		for (int j = 0; j < test_char[i]; j++)
	//		{
	//			fprintf(q_out, "%c", buffer_3[j]);
	//		}
	//		fprintf(q_out, "\n");
	//	}
	//	else
	//	{
	//		cout << "End of file. " << endl;
	//		for (int j = 0; j < length_1; j++)  // length-1：最后一位不是回车'\n'
	//		{
	//			fprintf(q_out, "%c", buffer_1[j]);
	//		}
	//		fprintf(q_out, "\n");
	//		for (int j = 0; j < test_char[i]; j++)
	//		{
	//			fprintf(q_out, "%c", buffer_2[j]);
	//		}
	//		fprintf(q_out, "\n");
	//		fprintf(q_out, "+");
	//		fprintf(q_out, "\n");
	//		for (int j = 0; j < test_char[i]; j++)
	//		{
	//			fprintf(q_out, "%c", buffer_3[j]);
	//		}
	//	}
	//	i++;
	//}


	delete[] acm_ny;
	delete[] quick;
	freearay(rows, input_hist);
	freearay(rows, at_cg);
	//fclose(qual_out);
	//fclose(q_out);
	fclose(fp);
}
int main(int argc, char* argv[])
{
	string file_str = argv[1];
	de_code(file_str);
}
