
// GREY_WOLF_OPTIMIZER.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

/*
Step1: Randomly initialize Grey wolf population of N particles Xi ( i=1, 2, �, n)


Step2: Calculate the fitness value of each individuals
           sort grey wolf population based on fitness values
           alpha_wolf = wolf with least fitness value
           beta_wolf  = wolf with second least fitness value
           gamma_wolf = wolf with third least fitness value
********************************************************************************************************

Step 3: For Iter in range(max_iter):  # loop max_iter times

            calculate the value of a
                a = 2*(1 - Iter/max_iter)

            For i in range(N):  # for each wolf

               a. Compute the value of A1, A2, A3 and C1, C2, C3
                     A1 = a*(2*r1 -1), A2 = a*(2*r2 -1), A3 = a*(2*r3 -1)
                     C1 = 2*r1, C2 = 2*r2, C3 = 2*r3

               b. Computer X1, X2, X3
                       X1 = alpha_wolf.position -
                             A1*abs(C1*alpha_wolf_position - ith_wolf.position)

                       X2 = beta_wolf.position -
                             A2*abs(C2*beta_wolf_position - ith_wolf.position)

                       X3 = gamma_wolf.position -
                             A3*abs(C3*gamma_wolf_position - ith_wolf.position)

               c. Compute new solution and it's fitness
                       Xnew = (X1 + X2 + X3) / 3
                       fnew = fitness( Xnew)

               d. Update the ith_wolf greedily
                     if( fnew < ith_wolf.fitness)

                         ith_wolf.position = Xnew
                         ith_wolf.fitness = fnew
             End-for

             # compute new alpha, beta and gamma
                   sort grey wolf population based on fitness values
                   alpha_wolf = wolf with least fitness value
                   beta_wolf  = wolf with second least fitness value
                   gamma_wolf = wolf with third least fitness value
         End -for



*/

/*
 * USERS GUIDE
 *
 * 1. EDIT THE BOUNDS -- (IN THE MAIN FUNCTION)
 * 2. EDIT THE OBJECTIVE FUNCTION -- (IN THE MAIN FUNCTION AND THE functions FUNCTION)
 * 3. EDIT THE POPULATION
 * 4. EDIT THE MAXIMUM NUMBER OF ITERATIONS -- TO MAYBE 100 OR MORE
 * */


#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <map>
#define MAX 10000000
#define MIN -1000

using namespace std;

const int population{ 5 }; //usually depicts number of wolves...
const int maxiter{ 100 }; //refers to the number of iterations
//The size of the array put down here is determined by how many factors you have
double functions(vector<double> wolf) {
    //Here you put your objective function...
//Step 4: Return best wolf in the population

//    double val;
//    double d = wolf.at(0);
//    double B = wolf.at(1);
//    double FCAR = wolf.at(2);
//    double t = wolf.at(3);
//    double FR = wolf.at(4);
//    double RS = wolf.at(5);

//    Objective functions for Ugo on the thermal drilling problem
    //Regression Equation

//    AF = 1.3693 + 0.01313 d - 0.00372 β - 0.00210 FCAR - 0.3009 t + 0.000329 FR + 0.000024 RS
//    val = 1.3693 + 0.01313*d - 0.00372 * B - 0.00210*FCAR - 0.3009 * t + 0.000329*FR+ 0.000024*RS;
//   RF	=	1.806 - 0.0860 d + 0.0023 β + 0.00074 FCAR - 0.1707 t - 0.00309 FR - 0.000040 RS
//    val = 1.806 - 0.0860 *wolf.at(0) +0.0023*wolf.at(1) + 0.00074*wolf.at(2) - 0.1707*wolf.at(3) - 0.00309*wolf.at(4) - 0.000040*wolf.at(5);
//    DE = 1.505 - 0.0730 d - 0.0051 β + 0.00025 FCAR - 0.090 t - 0.00402 FR + 0.000116 RS
//    val = 1.505 - 0.0730*wolf.at(0) - 0.0051*wolf.at(1) + 0.00025*wolf.at(2) - 0.090*wolf.at(3) - 0.00402*wolf.at(4) + 0.000116*wolf.at(5);
//    RE = 0.978 - 0.1405*d + 0.01184*B + 0.00616*FCAR - 0.3735*t + 0.00116*FR - 0.000053*RS;
//    val = 0.978 - 0.1405*wolf.at(0) + 0.01184*wolf.at(1)  + 0.00616*wolf.at(2)  - 0.3735*wolf.at(3)  + 0.00116*wolf.at(4)  - 0.000053*wolf.at(5);
    //BL =  -0.2265 + 0.08978 d - 0.000458 β - 0.001613 FCAR + 0.21064 t - 0.001149 FR + 0.000004 RS
//    val =  -0.2265 + 0.08978*d - 0.000458*B - 0.001613*FCAR + 0.21064*t - 0.001149*FR + 0.000004*RS;
//    val	=	1.292 + 0.0103*d - 0.00376*B - 0.00084*FCAR - 0.4090*t - 0.000475*FR + 0.000097*RS;
    double L = wolf.at(0);
    double S = wolf.at(1);
    double D = wolf.at(2);

    //OBJECTIVE FUNCTION FOR THE SPECIFIC WEAR RATE
//    double val = 0.000363 - 0.000008*L - 0.000160*S - 0.000000*D + 0.000000*L*L + 0.000035*S*S + 0.000000*D*D + 0.000002*L*S + 0.000000*L*D + 0.000000*S*D;
//    double val = 0.000228-(0.000002*L) - (0.000043*S) - (0.00000001*D);
      double val=0.3345 - 0.000386*L + 0.0574*S - 0.000013*D + 0.00000*L*L - 0.02100*S*S + 0.000000*D*D - 0.000928*L*S - 0.000000*L*D - 0.000000*S*D;
    return val;
}

vector<vector<double>> randomInit(vector<vector<double>> bounds) {
    //recall that the formular is x = L + r(U - L)
    //the aim is to return a vector with 5 rows and 7 columns, the last colunm being the f(all the wolves)
    //the vector bounds has the lower and upper boundaries...
    vector<vector<double>> init;
    vector<double>randoms{}; //Vector to store the random numbers for the first wolf

    double x;
    srand((unsigned)time(NULL)); //Seeding the random numbers so that they vary always

    for (int k{ 0 }; k < population; k++) {
        vector<double>init_part{};
        double random{};



        for (int i{ 0 }; i < bounds.size(); i++) {
            //generate the random numbers here
            random = (float)rand() / RAND_MAX; //generating random numbers between 0 and 1

            //Obtaining the random numbers only for the first wolf.
            if(k==0){
              randoms.push_back(random);
            }

            x = bounds.at(i).at(0) + random * (bounds.at(i).at(1) - bounds.at(i).at(0));
            init_part.push_back(x);
        }
        double fofx = functions(init_part);
        init_part.push_back(fofx); // you add the objective functions result of all the parameters into the vector
        init.push_back(init_part); // you push back one wolf into the initial eqn
    }
    cout<<"The random numbers generated for the first wolf is: \n";

    //displaying the random numbers for the first wolf.
    for(auto c:randoms){
        cout<<c<<"   ";
    }
    cout<<endl;

    return init;
}

void printVec(vector<vector<double>> wolf) {
    int wolfSize = wolf.at(0).size();
    for (int i{ 0 }; i < wolf.size(); i++) {
        for (int j{ 0 }; j < wolfSize; j++) {
            if (j == (wolfSize - 1)) {
                cout << "  ";
            }
            cout << setw(9) << left << wolf.at(i).at(j) << " ";
        }
        cout << endl;
    }
}

void printVecSingle(vector<double> wolf) {
    for (int j{ 0 }; j < wolf.size(); j++) {
        cout << setw(9) << left << wolf.at(j) << " ";
    }
    cout << endl;
}

vector<vector<double>> alphaDerive(vector<vector<double>> wolves) {
    //what you are comparing here is the last element of all the vectors in this vector of vectors.
    vector<double> alpha{ MAX };
    vector<double> beta{ MAX };
    vector<double> gamma{ MAX };
    vector < vector<double>> best{};

    int alpha_size = wolves.at(0).size();

    //filling the vectors with MAX< so that there is something to compare against ...
    for (int j{ 0 }; j < alpha_size; j++) {
        alpha.push_back(MAX);
        beta.push_back(MAX);
        gamma.push_back(MAX);
    }




    //this means that the elements you are comparing are in the wolves.at(i).at(alpha_size-1)
    for (int i{ 0 }; i < wolves.size(); i++) {
        if (wolves.at(i).at(alpha_size - 1) < alpha.at(alpha_size - 1)) {
            gamma = beta;
            beta = alpha;
            alpha = wolves.at(i);
        }
        else if (wolves.at(i).at(alpha_size - 1) < beta.at(alpha_size - 1)) {
            gamma = beta;
            beta = wolves.at(i);
        }
        else if (wolves.at(i).at(alpha_size - 1) < gamma.at(alpha_size - 1)) {
            gamma = wolves.at(i);
        }
    }

    //after the alpha, beta and gamma vectors must have been obtained
    best.push_back(alpha);
    best.push_back(beta);
    best.push_back(gamma);

    return best;
}

vector<vector<double>> alphaDerivemax(vector<vector<double>> wolves){
    //what you are comparing here is the last element of all the vectors in this vector of vectors.
        vector<double> alpha{ MIN };
        vector<double> beta{ MIN };
        vector<double> gamma{ MIN };
        vector < vector<double>> best{};

        int alpha_size = wolves.at(0).size();

        //filling the vectors with MAX< so that there is something to compare against ...
        for (int j{ 0 }; j < alpha_size; j++) {
            alpha.push_back(MIN);	//filling with zeros.............
            beta.push_back(MIN);
            gamma.push_back(MIN);
        }




        //this means that the elements you are comparing are in the wolves.at(i).at(alpha_size-1)

        for (int i{ 0 }; i < wolves.size(); i++) {
            if (wolves.at(i).at(alpha_size - 1) > alpha.at(alpha_size - 1)) {
                gamma = beta;
                beta = alpha;
                alpha = wolves.at(i);
            }
            else if (wolves.at(i).at(alpha_size - 1) > beta.at(alpha_size - 1)) {
                gamma = beta;
                beta = wolves.at(i);
            }
            else if (wolves.at(i).at(alpha_size - 1) > gamma.at(alpha_size - 1)) {
                gamma = wolves.at(i);
            }
        }

        //after the alpha, beta and gamma vectors must have been obtained
        best.push_back(alpha);
        best.push_back(beta);
        best.push_back(gamma);

        return best;
}

void graph(map<int, double> param) {
    cout << "THE BEST VALUES ARE AS FOLLOWS: \n";
    for (auto a : param) {
        cout << "Iteration " << a.first << " : " << a.second << endl;
    }
}

//make a function that runs the iterations. lets call it iterations hub for now.....
void iterations(vector<vector<double>> wolves, vector<vector<double>>best, vector<vector<double>>bounds) {
    /*
    Step 3: For Iter in range(max_iter):  # loop max_iter times
    */

    map<int, double> graph_param{};


    srand((unsigned)time(NULL));
    for (int i{ 1 }; i <= maxiter; i++) {
        int c_iterator = i;
        cout << "\t\t\n\nITERATION " << i << "\n";

        /*
        calculate the value of a
                a = 2*(1 - Iter/max_iter)

            For i in range(N):  # for each wolf

               a. Compute the value of A1, A2, A3 and C1, C2, C3
                     A1 = a*(2*r1 -1), A2 = a*(2*r2 -1), A3 = a*(2*r3 -1)
                     C1 = 2*r1, C2 = 2*r2, C3 = 2*r3
        */
        double a = 2 * (1 - i / maxiter); //a is redefined for each iteration...
        double r1 = (float)rand() / RAND_MAX;
        double r2 = (float)rand() / RAND_MAX;
        double r3 = (float)rand() / RAND_MAX;
        cout<<"r1: "<<r1<<endl;

        double A1, A2, A3;
        A1 = a * (2 * r1 - 1);
        A2 = a * (2 * r2 - 1);
        A3 = a * (2 * r3 - 1);

        r1 = (float)rand() / RAND_MAX;
        r2 = (float)rand() / RAND_MAX;
        r3 = (float)rand() / RAND_MAX;

        double C1, C2, C3;
        C1 = 2 * r1;
        C2 = 2 * r2;
        C3 = 2 * r3;

        cout<<"r1 for c1: "<<r1<<endl;
        /*

               b. Computer X1, X2, X3
                       X1 = alpha_wolf.position -
                             A1*abs(C1*alpha_wolf_position - ith_wolf.position)

                       X2 = beta_wolf.position -
                             A2*abs(C2*beta_wolf_position - ith_wolf.position)

                       X3 = gamma_wolf.position -
                             A3*abs(C3*gamma_wolf_position - ith_wolf.position)

                             c. Compute new solution and it's fitness
                       Xnew = (X1 + X2 + X3) / 3
                       fnew = fitness( Xnew)

               d. Update the ith_wolf greedily
                     if( fnew < ith_wolf.fitness)

                         ith_wolf.position = Xnew
                         ith_wolf.fitness = fnew
             End-for

             # compute new alpha, beta and gamma
                   sort grey wolf population based on fitness values
                   alpha_wolf = wolf with least fitness value
                   beta_wolf  = wolf with second least fitness value
                   gamma_wolf = wolf with third least fitness value
        */
        //loop through all the wolves here... finding xnew for each and carrying out the greedy selection
        for (int i{ 1 }; i <= population; i++) {
            cout << "\t\t\t\n\nWOLF " << i << "\n";
            //defining the alpha, beta and gamma wolves
            vector<double> alpha = best.at(0);
            alpha.pop_back(); //removing the last element which is the f (x) value...
            vector<double> beta = best.at(1);
            beta.pop_back();
            vector<double> gamma = best.at(2);
            gamma.pop_back();

            //FIND X1, X2, X3, XNEW, GREEDY SELECTION AND ON TO THE SECOND WOLF
            //
            //TO FIND X1
            vector<double> C1Xa;
            for (auto c : alpha) {
                c *= C1;
                C1Xa.push_back(c);
            }

            vector<double> Xt = wolves.at(i - 1);   //referring to the wolf being dealt with currently...
            Xt.pop_back();
            double previous = functions(Xt);


            vector<double> Dalpha;
            for (int i{ 0 }; i < Xt.size(); i++) {
                double val;
                val = abs(C1Xa.at(i) - Xt.at(i));
                Dalpha.push_back(val);
            }
            vector<double> X1;
            for (int j{ 0 }; j < alpha.size(); j++) {
                double vals = alpha.at(j) - (A1 * Dalpha.at(j));
                X1.push_back(vals);
            }

            //TO FIND X2
            vector<double> C2Xb;
            for (auto c : beta) {
                c *= C2;
                C2Xb.push_back(c);
            }

            //vector<double> Xt = wolves.at(i);   //referring to the wolf being dealt with currently...
            //Xt.pop_back();


            vector<double> Dbeta;
            for (int i{ 0 }; i < Xt.size(); i++) {
                double val;
                val = abs(C2Xb.at(i) - Xt.at(i));
                Dbeta.push_back(val);
            }
            vector<double> X2;
            for (int j{ 0 }; j < beta.size(); j++) {
                double vals = beta.at(j) - (A2 * Dbeta.at(j));
                X2.push_back(vals);
            }

            //TO FIND X3
            vector<double> C3Xg;
            for (auto c : gamma) {
                c *= C3;
                C3Xg.push_back(c);
            }

            //vector<double> Xt = wolves.at(i);   //referring to the wolf being dealt with currently...
            //Xt.pop_back();


            vector<double> Dgamma;
            for (int i{ 0 }; i < Xt.size(); i++) {
                double val;
                val = abs(C3Xg.at(i) - Xt.at(i));
                Dgamma.push_back(val);
            }
            vector<double> X3;
            for (int j{ 0 }; j < gamma.size(); j++) {
                double vals = gamma.at(j) - (A3 * Dgamma.at(j));
                X3.push_back(vals);
            }

            //To find Xnew
            vector<double> Xnew;
            for (int i{ 0 }; i < X1.size(); i++) {
                double valc;
                valc = (X1.at(i) + X2.at(i) + X3.at(i)) / 3;
                Xnew.push_back(valc);
            }

            cout << "\n\nX1: ";
            printVecSingle(X1);
            cout << "X2: ";
            printVecSingle(X2);
            cout << "X3: ";
            printVecSingle(X3);
            cout << "XN: ";
            printVecSingle(Xnew);

            //check if Xnew falls within the bounds,
            //If it does not, do not allow the greedy selection

            for (int i{ 0 }; i < Xnew.size(); i++) {
                for (int j{ 0 }; j <= 1; j++) {
                    if (Xnew.at(i) >= bounds.at(i).at(0) && Xnew.at(i) <= bounds.at(i).at(1)) {

                    }
                    else {
                        if (Xnew.at(i) < bounds.at(i).at(0)) {
                            Xnew.at(i) = bounds.at(i).at(0);
                        }
                        else if (Xnew.at(i) > bounds.at(i).at(1)) {
                            Xnew.at(i) = bounds.at(i).at(1);
                        }
                    }
                }
            }


            cout << "\n\n\t\tCarrying out Greedy Selection...\n";
            double latest = functions(Xnew);
            cout << "Previous Value: " << previous << endl;
            cout << "Xnew: " << latest << endl;
            if (previous > latest) { //greedy selection agreees for the Swap...
                cout << "Wolf will be updated\n";
                Xnew.push_back(latest);
                wolves.at(i - 1) = Xnew;
            }

            printVec(wolves);
        }

        //alpha, beta and gamma have to change with respect to your new wolves...
        best = alphaDerive(wolves);
        cout << "\n\tDISPLAYING NEW VALUES FOR ALPHA, BETA AND GAMMA\n\n";
        printVec(best);
        cout << endl;

        //Populate the map of the parameters to be printed here...
        int alpha_size = best.at(0).size();
        double alpha_value = best.at(0).at(alpha_size - 1);
        graph_param.insert(pair<int, double>(c_iterator, alpha_value));
    }
    graph(graph_param);
}

void iterationsmax(vector<vector<double>> wolves, vector<vector<double>>best, vector<vector<double>>bounds){
    /*
    Step 3: For Iter in range(max_iter):  # loop max_iter times
    */

    map<int, double> graph_param{};


    srand((unsigned)time(NULL));
    for (int i{ 1 }; i <= maxiter; i++) {
        int c_iterator = i;
        cout << "\t\t\n\nITERATION " << i << "\n";

        /*
        calculate the value of a
                a = 2*(1 - Iter/max_iter)

            For i in range(N):  # for each wolf

               a. Compute the value of A1, A2, A3 and C1, C2, C3
                     A1 = a*(2*r1 -1), A2 = a*(2*r2 -1), A3 = a*(2*r3 -1)
                     C1 = 2*r1, C2 = 2*r2, C3 = 2*r3
        */
        double a = 2 * (1 - i / maxiter); //a is redefined for each iteration...
        double r1 = (float)rand() / RAND_MAX;
        double r2 = (float)rand() / RAND_MAX;
        double r3 = (float)rand() / RAND_MAX;

        double A1, A2, A3;
        A1 = a * (2 * r1 - 1);
        A2 = a * (2 * r2 - 1);
        A3 = a * (2 * r3 - 1);

        double C1, C2, C3;
        C1 = 2 * r1;
        C2 = 2 * r2;
        C3 = 2 * r3;
        /*

               b. Computer X1, X2, X3
                       X1 = alpha_wolf.position -
                             A1*abs(C1*alpha_wolf_position - ith_wolf.position)

                       X2 = beta_wolf.position -
                             A2*abs(C2*beta_wolf_position - ith_wolf.position)

                       X3 = gamma_wolf.position -
                             A3*abs(C3*gamma_wolf_position - ith_wolf.position)

                             c. Compute new solution and it's fitness
                       Xnew = (X1 + X2 + X3) / 3
                       fnew = fitness( Xnew)

               d. Update the ith_wolf greedily
                     if( fnew < ith_wolf.fitness)

                         ith_wolf.position = Xnew
                         ith_wolf.fitness = fnew
             End-for

             # compute new alpha, beta and gamma
                   sort grey wolf population based on fitness values
                   alpha_wolf = wolf with least fitness value
                   beta_wolf  = wolf with second least fitness value
                   gamma_wolf = wolf with third least fitness value
        */
        //loop through all the wolves here... finding xnew for each and carrying out the greedy selection
        for (int i{ 1 }; i <= population; i++) {
            cout << "\t\t\t\n\nWOLF " << i << "\n";
            //defining the alpha, beta and gamma wolves
            vector<double> alpha = best.at(0);
            alpha.pop_back(); //removing the last element which is the f (x) value...
            vector<double> beta = best.at(1);
            beta.pop_back();
            vector<double> gamma = best.at(2);
            gamma.pop_back();

            //FIND X1, X2, X3, XNEW, GREEDY SELECTION AND ON TO THE SECOND WOLF
            //
            //TO FIND X1
            vector<double> C1Xa;
            for (auto c : alpha) {
                c *= C1;
                C1Xa.push_back(c);
            }

            vector<double> Xt = wolves.at(i - 1);   //referring to the wolf being dealt with currently...
            Xt.pop_back();
            double previous = functions(Xt);


            vector<double> Dalpha;
            for (int i{ 0 }; i < Xt.size(); i++) {
                double val;
                val = abs(C1Xa.at(i) - Xt.at(i));
                Dalpha.push_back(val);
            }
            vector<double> X1;
            for (int j{ 0 }; j < alpha.size(); j++) {
                double vals = alpha.at(j) - (A1 * Dalpha.at(j));
                X1.push_back(vals);
            }

            //TO FIND X2
            vector<double> C2Xb;
            for (auto c : beta) {
                c *= C2;
                C2Xb.push_back(c);
            }

            //vector<double> Xt = wolves.at(i);   //referring to the wolf being dealt with currently...
            //Xt.pop_back();


            vector<double> Dbeta;
            for (int i{ 0 }; i < Xt.size(); i++) {
                double val;
                val = abs(C2Xb.at(i) - Xt.at(i));
                Dbeta.push_back(val);
            }
            vector<double> X2;
            for (int j{ 0 }; j < beta.size(); j++) {
                double vals = beta.at(j) - (A2 * Dbeta.at(j));
                X2.push_back(vals);
            }

            //TO FIND X3
            vector<double> C3Xg;
            for (auto c : gamma) {
                c *= C3;
                C3Xg.push_back(c);
            }

            //vector<double> Xt = wolves.at(i);   //referring to the wolf being dealt with currently...
            //Xt.pop_back();


            vector<double> Dgamma;
            for (int i{ 0 }; i < Xt.size(); i++) {
                double val;
                val = abs(C3Xg.at(i) - Xt.at(i));
                Dgamma.push_back(val);
            }
            vector<double> X3;
            for (int j{ 0 }; j < gamma.size(); j++) {
                double vals = gamma.at(j) - (A3 * Dgamma.at(j));
                X3.push_back(vals);
            }

            //To find Xnew
            vector<double> Xnew;
            for (int i{ 0 }; i < X1.size(); i++) {
                double valc;
                valc = (X1.at(i) + X2.at(i) + X3.at(i)) / 3;
                Xnew.push_back(valc);
            }

            cout << "\n\nX1: ";
            printVecSingle(X1);
            cout << "X2: ";
            printVecSingle(X2);
            cout << "X3: ";
            printVecSingle(X3);
            cout << "XN: ";
            printVecSingle(Xnew);

            //check if Xnew falls within the bounds,
            //If it does not, do not allow the greedy selection

            for (int i{ 0 }; i < Xnew.size(); i++) {
                for (int j{ 0 }; j <= 1; j++) {
                    if (Xnew.at(i) >= bounds.at(i).at(0) && Xnew.at(i) <= bounds.at(i).at(1)) {

                    }
                    else {
                        if (Xnew.at(i) < bounds.at(i).at(0)) {
                            Xnew.at(i) = bounds.at(i).at(0);
                        }
                        else if (Xnew.at(i) > bounds.at(i).at(1)) {
                            Xnew.at(i) = bounds.at(i).at(1);
                        }
                    }
                }
            }


            cout << "\n\n\t\tCarrying out Greedy Selection...\n";
            double latest = functions(Xnew);
            cout << "Previous Value: " << previous << endl;
            cout << "Xnew: " << latest << endl;
            if (previous < latest) { //greedy selection agreees for the Swap...
                cout << "Wolf will be updated\n";
                Xnew.push_back(latest);
                wolves.at(i - 1) = Xnew;
            }

            printVec(wolves);
        }

        //alpha, beta and gamma have to change with respect to your new wolves...
        best = alphaDerivemax(wolves);
        cout << "\n\tDISPLAYING NEW VALUES FOR ALPHA, BETA AND GAMMA\n\n";
        printVec(best);
        cout << endl;

        //Populate the map of the parameters to be printed here...
        int alpha_size = best.at(0).size();
        double alpha_value = best.at(0).at(alpha_size - 1);
        graph_param.insert(pair<int, double>(c_iterator, alpha_value));
    }
    graph(graph_param);
}

int main(int argc, char *argv[])
{
   
    //Defining and initalising variables here sir...
      vector<double> alpha{};
      vector<double> beta{};
      vector<double> gamma{};
      vector<vector<double>>bestValues{};

      /*vector<vector<double>> bounds{ {1,4},{30,60},{10,22},{0.25,1} };*/ //my vector of boundaries. Each wolf has a boundary...

      //NUMBER 1 STEP
//      vector<vector<double>>bounds { {5.4, 9.2}, { 30,60 }, { 50,100 }, { 1,3 }, { 60,140 }, { 1500,3500 }}; //Boundaries for Ugo on Thermal Drilling problem
       vector<vector<double>>bounds { {10, 40}, { 0.42,1.68 }, { 750,3000 }};
      vector<vector<double>> randomInits = randomInit(bounds); //random initial values are generated here

      //Input your objective function here
//      cout<<"OBJECTIVE FUNCTION: BL = -0.2265 + 0.08978 d - 0.000458 β - 0.001613 FCAR + 0.21064 t - 0.001149 FR + 0.000004 RS\n";

      //ENTER YOUR OBJECTIVE FUNCTION
      cout<<"OBJECTIVE FUNCTION:μ	=	0.3345 - 0.000386 L + 0.0574 S - 0.000013 D + 0.000003 L*L - 0.02100 S*S + 0.000000 D*D - 0.000928 L*S - 0.000000 L*D - 0.000000 S*D\n";
      //Choose 1 or 2 to either minimize of maximize
      cout<<"Would you like to Maximize or minimize?\n1 - MAXIMIZE\n2 - MINIMIZE\nCHOICE: ";
      int choice{};
      cin>>choice;

      cout << "\tRandomly initialized Grey wolf population of " << population << " wolves\n\n";
      printVec(randomInits);


      cout << "\n\nCarrying out the Minimization Process\n\n";
      cout << "Step 1 --- Identifying the alpha, beta and gamma wolves\n\n";

    if (choice == 1){
       //Obtaining the Best values while maximizing...
        //Display the ALpha, Beta and Gamma as vectors in the main, and derive them from a function
        bestValues = alphaDerivemax(randomInits); //bestValues as we know should have 3 vectors... alpha, beta and gamma in that order
        printVec(bestValues);
        iterationsmax(randomInits, bestValues, bounds);
    }

    else{
        //Obtaining Best Values while minimizing
        //Display the ALpha, Beta and Gamma as vectors in the main, and derive them from a function
        bestValues = alphaDerive(randomInits); //bestValues as we know should have 3 vectors... alpha, beta and gamma in that order
        printVec(bestValues);

        //Construct a function that will carry out the iteration
        //the function will need the following as inputs:
        /*
        `````The vector that has the alpha beta and gamma wolves
             The vector with all the wolves at the moment
    */
        iterations(randomInits, bestValues, bounds);
    }



    return 0;
}
