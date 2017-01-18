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
		float t;	//σbk^2
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

	//约束目标
	for(int i = 0; i < 1; i++)
		BP += A[0][i];	//用A来定义了BP,以后BP用A代替
	model.add(IloMinimize(environment, BP));
	//约束1 Δσ>=σbk-σ~bk
	for(int i = 0; i < k; i++){
		IloExpr constraint(environment);
			constraint +=  S[0][i] - sigma_2_bk[i]- A[0][0];
		model.add(constraint <= 0);
	}
	//约束2 -Δσ<=σbk-σ~bk
	for(int i = 0; i < k; i++){
		IloExpr constraint(environment);
			constraint += S[0][i] - sigma_2_bk[i] + A[0][0];
		model.add(constraint >= 0);
	}
	//约束3 Σσbk = Kσ+ε
	IloExpr constraint(environment);
	for(int i = 0; i < k; i++){	
			constraint += S[0][i];		
	}
	constraint -= k*sigma+epis;
	model.add(constraint == 0);

	IloCplex solver(model);
	solver.solve();
	environment.out()<<"Δσ values   = "<<solver.getObjValue()<<endl;
	fprintf(logFile, "Δσ values = %f\n",solver.getObjValue());

	for(int i = 0; i < k; i++){
		float t = solver.getValue(S[0][i]);
		t = t*t;	//转回平方形式
		printf("σ%d = %f\n",i,t);
		ret.push_back(t);
		fprintf(logFile, "σ%d = %f\n",i,t);
	}
	fclose(logFile);
	return ret;
}


void main(){

	FILE *data = fopen("inputData//data.txt", "r");
	find_sigma(data);

}