#include <iostream>
#include <iomanip>
#include <math.h>

//#include "SCmathlib.h"


using namespace std;
//const int size;

int sizeInput, sizeN, size;


int main(){
    int count =0;
    bool flag;
    cout << "Enter the size of the system : "; 
    cin >> sizeInput;

    int maxMatsize = 25;
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
    for(i=0;i<sizeInput;i++){

        F[i] =  1.0;
        mainDia[i] = 2.0;
        subDia[i] = -1.0;
        supDia[i] = -1.0;
    }
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


/// Cyclic Reduction Step 1

for(i=0;i<log2(size+1)-1;i++){
    for(j=pow(2,i+1)-1;j<size;j=j+pow(2,i+1)){

        offset = pow(2,i);
       // cout << offset << "   " << j;
        index1 = j - offset;
        index2 = j + offset;

        alpha = subDia[j]/mainDia[index1];
        gamma = supDia[j]/mainDia[index2];

       // cout<< alpha << "   "<< gamma ; 
       // cout << endl;
    

    // if (j == size - 1){
    //     mainDia[j] = mainDia[j] - supDia[index1]*alpha;
    //     F[j] = F[j] - F[index1] * alpha;
    //     subDia[j] = -subDia[index1] * alpha;
    //     supDia[j] = 0.0;
    // }
    // else{
    //     mainDia[j] = mainDia[j] - subDia[index1] * alpha - supDia[index2] * gamma;
    //     F[j] = F[j] - F[index1] * alpha - F[index2] * gamma;
    //     subDia[j] = -subDia[index1] * alpha;
    //     supDia[j] = -supDia[index2] * gamma;
    // }

        subDia[j] = -subDia[index1]*(alpha);
        mainDia[j] = mainDia[j] - supDia[index1]*alpha - subDia[index2]*gamma;
        supDia[j] = -supDia[index2]*(gamma);
        F[j] = F[j] - F[index1] * alpha - F[index2] * gamma;
    }

}

    int index = (size - 1)/2;
    x[index] =  F[index]/mainDia[index];


    for(i=log2(size+1)-2;i>=0;i--){
        for(j=pow(2,i+1)-1;j<size;j=j+pow(2,i+1)){
            offset = pow(2,i);
            index1 = j - offset;
            index2 = j + offset;

            if (j != index1) {
            if (index1 - offset < 0){
               
                x[index1] = (F[index1]- supDia[index1]*x[index1+offset])/mainDia[index1];
            }
            else{
                x[index1] = (F[index1] - subDia[index1]*x[index1-offset] - supDia[index1]*x[index1+offset])/mainDia[index1];
            }
            }
        
        if(j != index2){
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
    return 0;  // return 0 to the OS.
}
