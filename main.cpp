#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include <ilcplex/ilocplex.h>

using namespace std;

vector<float> find_sigma(FILE *data) {
	FILE *logFile = fopen("outputData//logFile.txt", "w");
	vector<float> sigma_2_bk,ret;
	float sigma,epis;
	int k;
	fscanf(data,"%f %f %d", &sigma,&epis,&k);
	for(int i = 0; i < k; i++){
		float t;	//��bk^2
        fscanf(data, "%f", &t);
		t = sqrt(t);
		sigma_2_bk.push_back(t);
	}
	IloEnv environment;
    IloModel model(environment);

	IloArray<IloNumVarArray> A(environment,1);
	for(int i = 0; i < 1; i++)
		A[i] = IloNumVarArray(environment,1,0,99999);
	IloArray<IloNumVarArray> S(environment,1);
	for(int i = 0; i < 1; i++)
		S[i] = IloNumVarArray(environment,k,0,99999);
	IloExpr BP(environment);

	//Լ��Ŀ��
	for(int i = 0; i < 1; i++)
		BP += A[0][i];	//��A��������BP,�Ժ�BP��A����
	model.add(IloMinimize(environment, BP));
	//Լ��1 ����>=��bk-��~bk
	for(int i = 0; i < k; i++){
		IloExpr constraint(environment);
			constraint +=  S[0][i] - sigma_2_bk[i]- A[0][0];
		model.add(constraint <= 0);
	}
	//Լ��2 -����<=��bk-��~bk
	for(int i = 0; i < k; i++){
		IloExpr constraint(environment);
			constraint += S[0][i] - sigma_2_bk[i] + A[0][0];
		model.add(constraint >= 0);
	}
	//Լ��3 ����bk = K��+��
	IloExpr constraint(environment);
	for(int i = 0; i < k; i++){	
			constraint += S[0][i];		
	}
	constraint -= k*sigma+epis;
	model.add(constraint == 0);

	IloCplex solver(model);
	solver.solve();
	environment.out()<<"���� values   = "<<solver.getObjValue()<<endl;
	fprintf(logFile, "���� values = %f\n",solver.getObjValue());

	for(int i = 0; i < k; i++){
		float t = solver.getValue(S[0][i]);
		t = t*t;	//ת��ƽ����ʽ
		printf("��%d = %f\n",i,t);
		ret.push_back(t);
		fprintf(logFile, "��%d = %f\n",i,t);
	}
	fclose(logFile);
	return ret;
}


void main(){

	FILE *data = fopen("inputData//data.txt", "r");
	find_sigma(data);

}