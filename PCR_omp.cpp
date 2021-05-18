#include <iostream>
#include <iomanip>
#include <math.h>
#include <omp.h>
#include <time.h>

//#include "SCmathlib.h"


using namespace std;
//const int size;

int sizeInput, sizeN, size;

int maxMatsize = 25;


int main(){

    int count =0;
    bool flag;
    cout << "Enter the size of the system : "; 
    cin >> sizeInput;

    

    double * p = new double[maxMatsize];

     
    for (int k = 0; k<maxMatsize; k++ ){
        
        p[k] = pow(2,k)-1;

        if(p[k] == sizeInput){
            size = p[k];
            flag =  true;
            
           // cout << sizeN;
           // cout << endl;
        }
        else {
            if (p[k]>= sizeInput ){
                count = count +1;
            }

        }

        
    }

    
   
   
    sizeN = p[maxMatsize-count];
    if (flag ==  false){
        size = sizeN;
    }
   cout << size<< "  "<< sizeInput;
   cout << endl;



    
    clock_t start = clock();

    int i,j;
    int index1,index2,offset;
    double alpha,gamma;

    double * x = new double[size];
    for(i=0;i<size;i++){
    x[i] = 0.0;
    }

    double * F = new double[size];
    double *mainDia = new double[size];
    double *subDia = new double[size];
    double *supDia = new double[size];

    if (flag == false){    
    //#pragma omp parallel for
    for(i=0;i<sizeInput;i++){

        F[i] =  1.0;
        mainDia[i] = 2.0;
        subDia[i] = -1.0;
        supDia[i] = -1.0;
    }
    
    //#pragma omp parallel for 
    for (i=sizeInput; i< size; i++){

        F[i] =  1.0;
        mainDia[i] = 1.0;
        subDia[i] = 0.0;
        supDia[i] = 0.0;

    }
    
    
    }

    else{



    for(i=0;i<size;i++){

    
        F[i] =  1.0;
        mainDia[i] = 2.0;
        subDia[i] = -1.0;
        supDia[i] = -1.0;
    }
    
    }

    subDia[0] = 0.0;
    supDia[sizeInput-1] = 0.0;

int logSize = log2(size+1)-1;


/// Cyclic Reduction Step 1
for(i=0;i<logSize;i++){
    int step = pow(2,i+1);
    #pragma omp parallel shared (subDia, supDia, mainDia,F, size) private(j,index1, index2, alpha, gamma)
    {
    #pragma omp for
    for(j=pow(2,i+1)-1;j<size;j=j+ step){

        

        //offset = pow(2,i);
        index1 = j - pow(2,i);
        index2 = j + pow(2,i);

        alpha = subDia[j]/mainDia[index1];
        gamma = supDia[j]/mainDia[index2];

        
      
        //#pragma omp atomic capture
        subDia[j] = -subDia[index1]*(alpha);
        mainDia[j] = mainDia[j] - supDia[index1]*alpha - subDia[index2]*gamma;
        supDia[j] = -supDia[index2]*(gamma);
        F[j] = F[j] - F[index1] * alpha - F[index2] * gamma;
    }

    }

}

    int index = (size - 1)/2;
    x[index] =  F[index]/mainDia[index];


    for(i=log2(size+1)-2;i>=0;i--){
        int step = pow(2,i+1);
       #pragma omp parallel shared(x,F,subDia, supDia, mainDia, size) private(j,index1, index2, alpha, gamma)

       {
        #pragma omp  for    
        for(j=pow(2,i+1)-1;j<size;j=j+ step){
            offset = pow(2,i);
            index1 = j - offset;
            index2 = j + offset;

           // printf("Executed by %d \n", omp_get_thread_num());
            if (index1 - offset < 0){
               
                x[index1] = (F[index1]- supDia[index1]*x[index1+offset])/mainDia[index1];
            }
            else{
                x[index1] = (F[index1] - subDia[index1]*x[index1-offset] - supDia[index1]*x[index1+offset])/mainDia[index1];
            }
        
            if(index2 + offset >= size ){
                x[index2] = (F[index2] - subDia[index2]*x[index2-offset])/mainDia[index2];
            }
            else{
               x[index2] = (F[index2] - subDia[index2]*x[index2-offset] - supDia[index2]*x[index2+offset])/mainDia[index2];
            }
        }        

        }
    

    }


    
  for(i=0;i<sizeInput;i++){
        cout << x[i] << endl;
}

// Stop measuring time and calculate the elapsed time
    clock_t end = clock();
    double elapsed = double(end - start)/CLOCKS_PER_SEC;
    
//    printf("Time measured: %.3f seconds.\n", elapsed);
	cout << "The Time measured is : " << elapsed << endl;
    return 0;  // return 0 to the OS.
}
