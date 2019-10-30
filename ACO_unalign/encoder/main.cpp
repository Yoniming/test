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

#define Test_Line_maxsize 2000
int total_len = 100000;
static int min_0 = 100, max_0 = 0;
vector<int> test_char;
char *test_quick = new char[Test_Line_maxsize];
char test_input_buf[Test_Line_maxsize];
char **input_hist;
char **at_cg;
int rows = 0;
int allLen=0;
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


int get_max(string file_str)
{
	int length;
	const string file_path = "./" + file_str;
	FILE * test_fp;
	if (!(test_fp = fopen(file_path.c_str(), "r")))
	{
		cout << "Open Failed!" << endl;
		return -1;
	}
	while (!feof(test_fp))
	{
		fgets(test_quick, Test_Line_maxsize, test_fp);
		fgets(test_quick, Test_Line_maxsize, test_fp);
		fgets(test_quick, Test_Line_maxsize, test_fp);
		fgets(test_input_buf, Test_Line_maxsize, test_fp);
		length = strlen(test_input_buf) - 1;
		test_char.push_back(strlen(test_input_buf) - 1);
		allLen++;
		for (int i = 0; i < length; i++)
		{
			if (test_input_buf[i] < min_0)
			{
				min_0 = test_input_buf[i];
			}
			if (test_input_buf[i] > max_0)
			{
				max_0 = test_input_buf[i];
			}
		}
	}
	fclose(test_fp);
}
int compress(string file_str)
{
	char * tmp;
	char * mid;
	int line_max = *max_element(test_char.begin(), test_char.end());
	input_hist = array_char(total_len, line_max+2);
	at_cg = array_char(total_len, line_max + 2);
	int num_symbol = max_0 - min_0 + 1;

	int model_num = num_symbol*num_symbol*num_symbol * 16;
	int adaptive_flag = 1;
	/*int cols = line_max;*/
	char table[256];

	//===============================1.2 Encoder definition=========================//
	ac_encoder ace_ny;
	ac_model line_model[10];
	ac_model * acm_ny = new ac_model[model_num];
	memset(table, 0, 256);
	for (int i = 0; i < 10; i++)
		ac_model_init(&line_model[i], line_max+1, NULL, adaptive_flag);

	for (int i = 0; i < model_num; i++)
	{
		ac_model_init(&acm_ny[i], num_symbol, NULL, adaptive_flag);
	}

	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;

	const string file_out = file_str + ".stream";
	ac_encoder_init(&ace_ny, file_out.c_str());
	fwrite(&line_max, sizeof(int), 1, ace_ny.fp);
	fwrite(&max_0, sizeof(int), 1, ace_ny.fp);
	fwrite(&min_0, sizeof(int), 1, ace_ny.fp);
	fwrite(&allLen, sizeof(int), 1, ace_ny.fp);


	const string file_path = "./" + file_str;
	FILE * fp;
	if (!(fp = fopen(file_path.c_str(), "r")))
	{
		cout << "Open Failed!" << endl;
		return -1;
	}
	int loop = 0;


	//===========是否编码每行长度，有些数据中在第一行会给出，所以可以不用压缩！！！
	for (int i = 0; i < allLen-1; i++)
	{
		ac_encode_symbol(&ace_ny, &line_model[0], test_char[i]);
	}
	//=============================

	while (!feof(fp))
	{
		for (int i = 0; i < total_len; i++)
		{
			fgets(test_quick, Test_Line_maxsize, fp);
			mid = fgets(at_cg[i], test_char[loop*total_len + i] + 2, fp);
			fgets(test_quick, Test_Line_maxsize, fp);
			tmp = fgets(input_hist[i], test_char[loop*total_len + i] + 2, fp);
			if (tmp == NULL)
			{
				break;
			}
			for (int j = 0; j < line_max + 2; j++)
			{
				input_hist[i][j] = input_hist[i][j] - min_0;
			}
			rows++;
		}


		//for (int m = 0; m < test_char[3]; m++)
		//{
		//	cout << int(input_hist[3][m]) << " ";
		//}


		//for (int n = 0; n < test_char[0]; n++)
		//{
		//	cout << int(input_hist[0][n]) << " ";
		//}
		//for (int i = 0; i < rows; i++)
		//{
		//	for (int j = 0; j < test_char[loop*total_len + i]; j++)
		//	{
		//		if (int(input_hist[i][j])<0)
		//		{
		//			cout << "the i is : " << i<<endl;
		//			cout << int(input_hist[i][j]) << " ";
		//		}
		//	}
		//}
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
					ac_encode_symbol(&ace_ny, &acm_ny[model_idx], input_hist[i][j]);
				}
			}
		rows = 0;
		loop++;
	}
	ac_encoder_done(&ace_ny);
	for (int i = 0; i < model_num; i++)
	{
		ac_model_done(&acm_ny[i]);
	}
	delete[] acm_ny;
	delete[] test_quick;
	freearay(rows, input_hist);
	freearay(rows, at_cg);
	fclose(fp);
}
int main(int argc, char* argv[])
{
	string file_str = argv[1];
	get_max(file_str);
	//vector<int>::iterator line_value;
	//for (line_value = test_char.begin(); line_value!=test_char.end(); line_value++)
	//{
	//cout << test_char[1600]<<" "<<test_char[1500] << " ";
	//}
	//cout << endl << " allLen: " << allLen << endl;
	//cout << *max_element(test_char.begin(), test_char.end()) << endl;
	//cout << min_0 << endl << max_0;
	//for (int i = 0; i < 1600; i++)
	//{
	//	cout << " test_char[" << i << "]" << "  " << test_char[i] << "   ";
	//}
	compress(file_str);
}
