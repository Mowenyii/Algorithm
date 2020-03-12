//这是一个非常简单的遗传算法源代码
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
using namespace std;

#define POPSIZE 100/* population size 种群大小*/ 
#define PXOVER 0.7/* probability of crossover 杂交概率*/ 
#define PMUTATION 0.07/* probability of mutation 变异概率*/ 
#define PI 3.14159265358979323846
#define MAXGENS 2000/* max. number of generations 最大代数*/

int generation; /* current generation no. 当前基因个数*/  
int cur_best; /* best individual 最优个体*/  

struct genotype
{
	double gene[2];/* a string of variables 变量*/ 
	double fitness;/* GT's fitness 基因的适应度*/
    double rfitness; /* relative fitness 比较适应度*/   
    double cfitness; /* cumulative fitness 积累适应度*/
};

struct genotype population[POPSIZE+1]; /* population 种群*/   
struct genotype newpopulation[POPSIZE+1]; /* new population; 新种群*/  /* replaces the old generation */  //取代旧的基因

/* Declaration of procedures used by this genetic algorithm */  
//以下是一些函数声明 
void initialize(); //种群基因结构体初始化
double randval(double,double);//随机数产生函数
void evaluate();//评价函数，可以由用户自定义，该函数取得每个基因的适应度
void keep_the_best();//保存每次遗传后的最佳基因
void elitist();//搜寻杰出个体函数：找出最好和最坏的个体。如果某代的最好个体比前一代的最好个体要坏，那么后者将会取代当前种群的最坏个体 
void select();//选择函数，轮盘赌
void crossover(double);//杂交函数：选择两个个体来杂交，这里用单点杂交
void Xover(int,int);//交叉 
void mutate(double);//变异函数：被该函数选中后会使得某一变量被一个随机的值所取代 

/////////////////////////////////////////////////////////////////
void initialize(){
	for(int j = 0; j < POPSIZE;j++){
		population[j].cfitness=0;//积累适应度
		population[j].fitness=0; //基因的适应度
		population[j].rfitness=0;//比较适应度
		population[j].gene[0]=randval(-3.0,12.1);
		population[j].gene[1]=randval(4.1,5.8);
	}	
}

double randval(double low,double high){
	double val;
	val = ((double)(rand()%1000)/1000.0)*(high-low)+low;
	return val;
}

//本题设两个基因x1,x2
//评价函数，根据题目要求f(x1,x2)=21.5+x1*sin(4*PI*x1)+x2*sin(20*PI*x2)
void evaluate(){
	for(int mem=0;mem<POPSIZE;mem++){
		double x1 = population[mem].gene[0];
		double x2 = population[mem].gene[1];
		population[mem].fitness=21.5+x1*sin(4*PI*x1)+x2*sin(20*PI*x2);
	}
}

//保存每次遗传后的最佳基因
void keep_the_best(){
	cur_best = 0;
	population[POPSIZE].fitness = population[0].fitness;
	population[POPSIZE].gene[0]=population[0].gene[0];
	population[POPSIZE].gene[1]=population[0].gene[1];
	for(int mem = 1;mem<POPSIZE;mem++){
		if(population[mem].fitness>population[POPSIZE].fitness){
			cur_best = mem;
			population[POPSIZE].fitness=population[mem].fitness;
		}
	}
	//一旦找到种群的最佳个体，就拷贝他的基因
	population[POPSIZE].gene[0]=population[cur_best].gene[0];
	population[POPSIZE].gene[1]=population[cur_best].gene[1];
}

//搜寻杰出个体函数：找出最好和最坏的个体。  
//如果某代的最好个体比前一代的最好个体要坏，那么后者将会取代当前种群的最坏个体 
void elitist(){
	double best,worst;//最好和最坏个体的适应度值
	int best_mem,worst_mem;//最好和最坏个体的 索引
	best = population[0].fitness;
	worst = best;
	for(int i=0;i<POPSIZE-1;i++){
		if(population[i].fitness > population[i+1].fitness){
			if(population[i].fitness >= best){
				best = population[i].fitness;
				best_mem = i;
			}
			if(population[i+1].fitness <= worst){
				worst = population[i+1].fitness;
				worst_mem = i + 1;
			}
		}
		else{
			if(population[i].fitness <= worst){
				worst = population[i].fitness;
				worst_mem = i;
			}
			if(population[i+1].fitness >= best){
				best = population[i+1].fitness;
				best_mem = i + 1;
			}
		}
	}
	//如果新种群中的最好个体比前一代的最好个体要强的话，那么就把新种群的最好个体拷贝出来。  
    //否则就用前一代的最好个体取代这次的最坏个体 
	if(best >= population[POPSIZE].fitness){
		population[POPSIZE].fitness=population[best_mem].fitness;
		population[POPSIZE].gene[0]=population[best_mem].gene[0];
		population[POPSIZE].gene[1]=population[best_mem].gene[1];
	}
	else{
		population[worst_mem].fitness=population[POPSIZE].fitness;
		population[worst_mem].gene[0]=population[POPSIZE].gene[0];
		population[worst_mem].gene[1]=population[POPSIZE].gene[1];
	}
}

//选择函数：用于最大化合并杰出模型的标准比例选择，保证最优秀的个体得以生存 
//轮盘赌选择法
void select(){
	double sum = 0.0;
	double p;
	int mem,i,j;
	/* find total fitness of the population */   
    //找出种群的适应度之和
	for(mem = 0;mem<POPSIZE;mem++){
		sum += population[mem].fitness;
	}
	/* calculate relative fitness */   
    //计算相对适应度
	for(mem = 0;mem<POPSIZE;mem++){
		population[mem].rfitness = population[mem].fitness/sum;
	}
	/* calculate cumulative fitness */   
    //计算累加适应度
	population[0].cfitness = population[0].rfitness;
	
	for(mem = 1;mem<POPSIZE;mem++){
		population[mem].cfitness = population[mem-1].cfitness+population[mem].rfitness;
	}
	/* finally select survivors using cumulative fitness. */   
    //用累加适应度作出选择
	for(i = 0;i<POPSIZE;i++){
		p = (double)(rand()%1000/1000.0);
		if(p < population[0].cfitness)
			newpopulation[i] = population[0];
		else{
			for(j = 0;j<POPSIZE;j++){
				if(p >= population[j].cfitness && p < population[j+1].cfitness)
					newpopulation[i] = population[j+1];
			}
		}
	}
	/* once a new population is created, copy it back */   
    //当一个新种群建立的时候，将其拷贝回去
	for(i=0;i<POPSIZE;i++){
		population[i]=newpopulation[i];
	}
}

//杂交函数：选择两个个体来杂交，这里用单点杂交 
void crossover(double pc){
	int mem,one;
	int first = 0;
	double x;
	for(mem = 0;mem < POPSIZE;mem++){
		x = (double)(rand()%1000/1000.0);
		if(x < pc){
			++first;
			if(first % 2==0)
				Xover(one,mem);
			else 
				one = mem;
		}
	}
}

//交叉   
void Xover(int one,int two){
	double temp;
	temp = population[one].gene[0];
	population[one].gene[0] = population[two].gene[0];
	population[two].gene[0] = temp;
}

//变异函数：被该函数选中后会使得某一变量被一个随机的值所取代 
void mutate(double pm){
	int i;
	double p1,p2;
	for(i=0;i<POPSIZE;i++){
		p1 = (double)(rand()%1000/1000.0);
		p2 = (double)(rand()%1000/1000.0);
		if(p1 < pm)
			population[i].gene[0]=randval(-3,12.1);
		if(p2 < pm)
			population[i].gene[1]=randval(4.1,5.8);
	}
}
/*
int main(){
	int i,j;
	for(i=1;i<10;i+=1){
		for(j=1;j<10;j+=1){
	initialize();
	evaluate();
	keep_the_best();
	generation=0;
	while(generation < MAXGENS){
		generation++;
		select();
		crossover(double(1.0*i/10));
		mutate(double(1.0*j/100));
		evaluate();
		elitist();
	}
	cout<<"0."<<i<<" 0.0"<<j<<" ";
	//cout<<setprecision(20)<<population[POPSIZE].gene[0]<<" "<<population[POPSIZE].gene[1]<<endl;
	cout<<setprecision(20)<<population[POPSIZE].fitness<<endl;
		}
	}
	return 0;
}
*/
int main()
{
	int i;
	
	for (int k = 1; k <= 9; k++)
	{
		for (int j = 1; j <= 9; j++)
		{
			double sum = 0.0;
			for (i = 1; i <= 30; i += 1)
			{
				initialize();
				evaluate(); //评价函数，可以由用户自定义，该函数取得每个基因的适应度
				keep_the_best();//保存每次遗传后的最佳基因
				generation = 0;
				while (generation < MAXGENS)
				{
					generation++;
					select();//选择函数：用于最大化合并杰出模型的标准比例选择，保证最优秀的个体得以生存		
					crossover(1.0 * k / 10);//杂交函数：选择两个个体来杂交，这里用单点杂交 
					//交配概率为 1.0 * k / 10	
					mutate(1.0 * j / 100);//变异函数：被该函数选中后会使得某一变量被一个随机的值所取代 
					//变异概率为1.0 * j / 10
					evaluate();//评价函数，可以由用户自定义，该函数取得每个基因的适应度
					elitist();//搜寻杰出个体函数：找出最好和最坏的个体。如果某代的最好个体比前一代的最好个体要坏，那么后者将会取代当前种群的最坏个体 
				}
				sum += population[POPSIZE].fitness;
				cout<<setprecision(20)<<population[POPSIZE].gene[0]<<" "<<population[POPSIZE].gene[1]<<endl;
				//cout<<setprecision(20)<<population[POPSIZE].fitness<<endl;
			}
			//cout<<setprecision(20)<<population[POPSIZE].gene[0]<<" "<<population[POPSIZE].gene[1]<<endl;
			//cout <<setprecision(3)<< 1.0 * k / 10 << " " << 1.0 * j / 100 << " " <<setprecision(20)<< sum / 30<<endl;
		}
	}
	return 0;
}