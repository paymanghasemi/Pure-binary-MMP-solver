/* 
 * File:   main.cpp
 * Author: c3156840
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ilcplex/ilocplex.h>
#include <sys/time.h>
#include <ctime>
#include <math.h>   
#include <vector>

using namespace std;

/*------------------------------------------------------------------------------
 
 -----------------------------Public Variables----------------------------------
 
 -----------------------------------------------------------------------------*/
#define Output_Precision 20
#define CPLEX_Relative_Gap 1.0e-6
#define Num_threads 1 
#define Time_Limit 3600 
#define Negative_infinity -1000000
#define Positive_infinity 1000000
#define epsilon 1.0e-6
#define Abs_Optimal_Tol 1.0e-6
#define Rel_Optimal_Tol 1.0e-6
#define denominator 1.0e-12


char* IP_file_name;
int N_Variables;
int N_Objectives;
int N_Constraints;
struct timeval startTime, endTime;
double totalTime(0);
double clock_Run_time(0);
double clock_start_time;
double run_start = 0;
double run_time = 0;
double T_limit = 0;

double* Y_star;
double Nodes_explored = 0;
double Nodes_proned = 0;
double Nodes_infeasible = 0;
double Global_LB;
double Global_UB;
double N_LP = 0;
double Check_RHS;
double Nodes_N_added = 0;
/*------------------------------------------------------------------------------
 
 ------------------------Declaring CPLEX information----------------------------
 
 -----------------------------------------------------------------------------*/
ILOSTLBEGIN
IloEnv env;
IloModel Model(env);
IloObjective Objective(env);
IloRangeArray Constraints(env);
IloSOS1Array sos1(env);
IloSOS2Array sos2(env);
IloNumVarArray Variable(env);
IloExprArray ObjF(env);
IloCplex Cplex(Model);
IloRangeArray Cuts(env);
IloRangeArray Tabulist(env);
IloRangeArray Bounds(env);
IloExpr Math(env);

IloModel Check_Model(env);
IloObjective Check_Objective(env);
IloCplex Check_Cplex(Check_Model);

/*------------------------------------------------------------------------------
 
 ------------------------------------Functions-----------------------------------
 
 -----------------------------------------------------------------------------*/


void Reading_LP_File_and_Generating_CPLEX_Variables() {
    Cplex.importModel(Model, IP_file_name, Objective, Variable, Constraints, sos1, sos2);
    sos1.clear();
    sos1.end();
    sos2.clear();
    sos2.end();
    Model.remove(Objective);
    N_Variables = Variable.getSize();
    ObjF.end();
    ObjF = IloExprArray(env, N_Objectives);
    Y_star = new double [N_Objectives];
    Model.remove(Constraints);
    for (int i = 0; i < N_Objectives; i++) {
        Y_star = 0;
        ObjF[i] = Constraints[0].getExpr();
        Constraints.remove(0);
    }
    N_Constraints = Constraints.getSize() - N_Objectives;
    Model.add(Constraints);
    Check_Model.add(Constraints);
    Math.clear();
    for (int i = 0; i < N_Objectives; i++) {
        Math += ObjF[i];
    }
    Objective = IloMaximize(env, Math);
    Model.add(Objective);

    Check_RHS = N_Objectives;
}

class Node {
public:

    double* Y;
    double Objective_value;
    double* LB_for_Obj;
    bool Infeasible;
    double UB;
    double LB;
    int Identifier;
    bool Do_Branch;

    Node() {
        LB_for_Obj = new double [N_Objectives];
        Y = new double [N_Objectives];
        for (int i = 0; i < N_Objectives; i++) {
            LB_for_Obj[i] = Negative_infinity;
        }
        Infeasible = 0;
    }

    void Reinitializing_The_Node(Node* Parent, int identifier) {
        for (int i = 0; i < N_Objectives; i++) {
            LB_for_Obj[i] = Parent->LB_for_Obj[i];
        }
        Identifier = identifier;
        LB_for_Obj[Identifier] = Parent->Y[Identifier];
        UB = Parent->UB;
    }

    bool Check_cutting_plane() {
        bool Check_Cutting_Plane_Value = 1;
        double Check_Value = 0;
        Math.clear();
        for (int i = 0; i < N_Objectives; i++) {
            Math += (ObjF[i] / Y[i]);
        }
        Check_Objective = IloMaximize(env, Math);
        Math.clear();

        Check_Model.add(Check_Objective);
        Check_Model.add(Tabulist);
        Check_Model.add(Cuts);
        Check_Model.add(Bounds);
        Check_Cplex.extract(Check_Model);
        Check_Cplex.setOut(env.getNullStream());
        Check_Cplex.setParam(IloCplex::Threads, Num_threads);
        Check_Cplex.setParam(IloCplex::TiLim, T_limit);
        Check_Cplex.setParam(IloCplex::EpGap, CPLEX_Relative_Gap);
        N_LP++;
        if (Check_Cplex.solve()) {
            Check_Value = Check_Cplex.getObjValue();
        }
        Check_Model.remove(Check_Objective);
        Check_Model.remove(Tabulist);
        Check_Model.remove(Cuts);
        Check_Model.remove(Bounds);
        Check_Model.remove(Check_Objective);

        if (Check_Value <= Check_RHS) {
            Check_Cutting_Plane_Value = 0;
        }


        return Check_Cutting_Plane_Value;
    }

    void Cutting_plane() {
        Do_Branch = Check_cutting_plane();
        Math.clear();
        for (int i = 0; i < N_Objectives; i++) {
            Math += (ObjF[i] / Y[i]);
        }
        Cuts.add(Math - Check_RHS >= 0);
        Math.clear();
    }

    void Node_Explore() {
        run_start = clock();
        Nodes_explored++;
        Bounds.clear();
        for (int i = 0; i < N_Objectives; i++) {
            Bounds.add(ObjF[i] >= LB_for_Obj[i] - epsilon);
        }
        Model.add(Tabulist);
        Model.add(Cuts);
        Model.add(Bounds);
        Cplex.extract(Model);
        Cplex.setOut(env.getNullStream());
        Cplex.setParam(IloCplex::Threads, Num_threads);
        Cplex.setParam(IloCplex::TiLim, T_limit);
        Cplex.setParam(IloCplex::EpGap, CPLEX_Relative_Gap);
        N_LP++;
        if (Cplex.solve()) {
            Objective_value = Cplex.getObjValue();
            LB = 1;
            for (int i = 0; i < N_Objectives; i++) {
                Y[i] = Cplex.getValue(ObjF[i]);
                LB *= Y[i];
            }
            UB = Objective_value / N_Objectives;
            UB = pow(UB, N_Objectives);

            Math.clear();
            for (int i = N_Objectives; i < N_Variables; i++) {
                if (Cplex.getValue(Variable[i]) >= 1 - epsilon) {
                    Math += (1 - Variable[i]);
                } else {
                    Math += Variable[i];
                }
            }
            Tabulist.add(Math >= 1);
            Math.clear();
            Do_Branch = 1;
        } else {
            Infeasible = 1;
            Nodes_infeasible++;
        }
        run_time = clock() - run_start;
        run_time = run_time / CLOCKS_PER_SEC;
        T_limit -= run_time;
        if (T_limit < 0) {
            T_limit = 0;
        }

        run_start = clock();
        if (Infeasible == 0 && LB > epsilon) {
            Cutting_plane();
        }
        Model.remove(Tabulist);
        Model.remove(Bounds);
        Model.remove(Cuts);
        Bounds.clear();

        run_time = clock() - run_start;
        run_time = run_time / CLOCKS_PER_SEC;
        T_limit -= run_time;
        if (T_limit < 0) {
            T_limit = 0;
        }
    }

    virtual ~Node() {
        delete [] LB_for_Obj;
        delete [] Y;
    }
};

vector <Node*>Tree_of_Nodes;

void Writing_Output() {
    ofstream OutputFile;
    OutputFile.open("Output.txt", ios::app);
    OutputFile << IP_file_name << " GLB= " << std::setprecision(10) << Global_LB << " GUB= " << Global_UB << " Gap= " << (Global_UB - Global_LB) / Global_UB << " #VAR= " << N_Variables - N_Objectives << " #Const= " << N_Constraints << " Time= " << (clock_Run_time / CLOCKS_PER_SEC) << " Open_nodes= " << Tree_of_Nodes.size() << " Explored_nodes= " << Nodes_explored << " Infeasible_nodes= " << Nodes_infeasible << " Pruned_nodes= " << Nodes_proned << " Nodes_not_added= " << Nodes_N_added*N_Objectives << " #IP= "<<N_LP<<endl;
    OutputFile.close();
}

void Add_The_Node_To_Tree(Node* New_Generated_Node) {
    run_start = clock();
    Node * New_Node[N_Objectives];

    //*********//
    for (int i = 0; i < N_Objectives; i++) {
        New_Node[i] = new Node;
        New_Node[i]->Reinitializing_The_Node(New_Generated_Node, i);
    }


    bool It_is_Added = 0;
    for (int i = 1; i < Tree_of_Nodes.size(); i++) {
        if (New_Generated_Node->UB > Tree_of_Nodes.at(i)->UB + epsilon) {
            for (int j = 0; j < N_Objectives; j++) {
                Tree_of_Nodes.insert(Tree_of_Nodes.begin() + i, New_Node[j]);
            }

            It_is_Added = 1;
            break;
        }
    }
    if (It_is_Added == 0) {
        for (int j = 0; j < N_Objectives; j++) {
            Tree_of_Nodes.push_back(New_Node[j]);
        }
    }
    run_time = clock() - run_start;
    run_time = run_time / CLOCKS_PER_SEC;
    T_limit -= run_time;
    if (T_limit < 0) {
        T_limit = 0;
    }
}

void Branch_and_bound() {
    while (T_limit > 0 && Tree_of_Nodes.size() > 0 && (Global_UB > Global_LB + Abs_Optimal_Tol) && (Global_UB - Global_LB) / (Global_UB + denominator) > Rel_Optimal_Tol) {

        Tree_of_Nodes.front()-> Node_Explore();
        if (Tree_of_Nodes.front()->Infeasible == 0 && Tree_of_Nodes.front()->UB > Global_LB) {
            if (Global_LB < Tree_of_Nodes.front()->LB) {
                Global_LB = Tree_of_Nodes.front()->LB;
            }
            if (Tree_of_Nodes.front()->Do_Branch == 1) {
                Add_The_Node_To_Tree(Tree_of_Nodes.front());
            } else {
                Nodes_N_added++;
            }
        } else if (Tree_of_Nodes.front()->Infeasible == 0 && Tree_of_Nodes.front()->UB <= Global_LB) {
                Nodes_proned++;
            }
        
        Tree_of_Nodes.front()->~Node();
        Tree_of_Nodes.erase(Tree_of_Nodes.begin());
        Global_UB = Tree_of_Nodes.front()->UB;

        if (T_limit > epsilon && (Global_UB <= Global_LB + Abs_Optimal_Tol || Tree_of_Nodes.size() == 0 || (Global_UB - Global_LB) / (Global_UB + denominator) <= Rel_Optimal_Tol)) {
            Global_UB = Global_LB;
        }
    }
}

int main(int argc, char *argv[]) {
    //---------------------------Preparation Phase------------------------------
    IP_file_name = argv[1];
    N_Objectives = atoi(argv[2]);
    T_limit = Time_Limit;
    Reading_LP_File_and_Generating_CPLEX_Variables();
    Global_LB = Negative_infinity;
    Global_UB = Positive_infinity;
    gettimeofday(&startTime, NULL);
    clock_start_time = clock();

    Node* Initial_Node = new Node;
    Tree_of_Nodes.push_back(Initial_Node);
    Branch_and_bound();

    gettimeofday(&endTime, NULL);
    clock_Run_time += (clock() - clock_start_time);
    totalTime += ((endTime.tv_sec - startTime.tv_sec) * 1000000L);
    totalTime += (endTime.tv_usec - startTime.tv_usec);
    Writing_Output();
    return 0;
}

