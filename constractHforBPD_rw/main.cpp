#include <iostream>
#include <fstream>
#include<cmath>
#include <itpp/itbase.h>
#include <chrono>
#include "func.hpp" 

//To run: make bike
//        ./bike.out

int main(){

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    int tempcof, rank, d, Ad;
    long coff, g_deg;
    long verbose = 0;
    //To change: l = n/2
    int l = 31;
    //
    itpp::bin tempcof2;
    itpp::GF2mat Aa(l,l), Bb(l, l), TAa(l,l), TBb(l,l), Hx(l,2*l), Hz(l,2*l);
    itpp::bvec vtemp(l), a(l), b(l);

    //construct polinomial a
    //To change: a(x) = ...
    a.zeros();
    a(0) = 1;
    a(1) = 1;
    a(3) = 1;
    //
    Aa = circulant_mat(a); 

    

    //construct polynomial b
    //To change: b(x) = ...
    b.zeros();
    b(0) = 1;
    b(4) = 1;
    b(17) = 1;
    //
    Bb = circulant_mat(b);

    //construct Hx
    for(int i = 0; i < l; i++){
        vtemp = Aa.get_col(i);
        Hx.set_col(i, vtemp);
    }
    for(int j = l; j < 2*l; j++){
        vtemp = Bb.get_col(j-l);
        Hx.set_col(j, vtemp);
    }

    //construct Hz
    TAa = Aa.transpose();
    TBb = Bb.transpose();
    for(int p = 0; p < l; p++){
        vtemp = TBb.get_col(p);
        Hz.set_col(p, vtemp);
    }
    for(int q = l; q < 2*l; q++){
        vtemp = TAa.get_col(q-l);
        Hz.set_col(q, vtemp);
    }

    std::cout << "for xingrui:" << "\n";
    std::cout << "Hx" << "\n";
    std::cout << 2*l << " " << l << "\n";

    for(int onei = 0; onei < l; onei++){
        for(int onej = 0; onej < 2*l; onej++){
            if(Hx(onei,onej) == 1){
                std::cout << onej+1 << " ";
            }
        }
        std::cout << "\n"; 
    }

    std::cout << "Hz" << "\n";
    std::cout << 2*l << " " << l << "\n";

    for(int onek = 0; onek < l; onek++){
        for(int onel = 0; onel < 2*l; onel++){
            if(Hz(onek,onel) == 1){
                std::cout << onel+1 << " ";
            }
        }
        std::cout << "\n"; 
    }

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
    std::cout << "Run-time " << time_span.count() << " seconds.\n";
    std::cout << "Run-time " << time_span.count()/60 << " minutes.\n";
    std::cout << "Run-time " << time_span.count()/3600 << " hours.\n";

    return 0;
}