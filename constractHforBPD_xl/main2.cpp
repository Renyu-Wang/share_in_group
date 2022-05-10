#include <iostream>
#include <fstream>
#include<cmath>
#include <itpp/itbase.h>
#include <chrono>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/GF2XFactoring.h>
#include "func.hpp"
#include "weilei_lib/weilei_lib.h" 

//To run: make bike
//        ./bike.out

int main(int argc, char ** argv){

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    int tempcof, rank, d, Ad;
    long coff, g_deg;
    long verbose = 0;
    //To change: l = n/2
    
    std::istringstream argv1( argv[1] );
    int l;
    if (argv1>>l){}
    else{std::cout<<"l need to be an int"<<std::endl;return 1;}

    std::istringstream argv2( argv[2] );
    int a_len;
    if (argv2>>a_len){}
    else{std::cout<<"a_len need to be an int"<<std::endl;return 1;}

    std::istringstream argv3( argv[3] );
    int b_len;
    if (argv3>>b_len){}
    else{std::cout<<"b_len need to be an int"<<std::endl;return 1;}

    if (argc!=a_len+b_len+4){std::cout<<"the number of input parameters is not right, you input"<<argc-1<<"parameters,  but a_len+b_len+3="<<a_len+b_len+3<<std::endl;return 1;}
    
    //
    itpp::bin tempcof2;
    NTL::GF2 c,z;
    itpp::GF2mat Aa(l,l), Bb(l, l), TAa(l,l), TBb(l,l), Hx(l,2*l), Hz(l,2*l);
    itpp::bvec vtemp(l), a(l), b(l);
    NTL::GF2X l_pol, a_pol, b_pol, g1, g2, m, c1_pol, c2_pol, x1, x2, ab_pol, c_pol, y;
    //NTL::vec_GF2 roots;
    NTL::vec_pair_GF2X_long l_factors, a_factors, b_factors;

    //construct polinomial l
    l_pol = NTL::GF2X(l,1) + 1;

    //construct polinomial a
    //To change: a(x) = ...
    a_pol += NTL::GF2X(1,1) + NTL::GF2X(3,1) + 1;
    a.zeros();
    // a(0) = 1;
    //a(1) = 1;
    // a(3) = 1;
    //
    for (int i=0;i<a_len;i++)
      {
	std::istringstream argvi( argv[i+4] );
	 if (argvi>>a(i)){}
	 else {std::cout<<"a(i) need to be binary"<<std::endl;return 1;}
      }
    Aa = circulant_mat(a); 

    

    //construct polynomial b
    //To change: b(x) = ...
    b_pol += NTL::GF2X(1,1) + NTL::GF2X(3,1) + 1;
    //b_pol =1 ;
    b.zeros();
      for (int i=0;i<b_len;i++)
      {
	std::istringstream argvi( argv[i+4+a_len] );
	 if (argvi>>b(i)){}
	 else {std::cout<<"b(i) need to be binary"<<std::endl;return 1;}
      }
    // b(0) = 1;
    // b(4) = 1;
    //b(17) = 1;
    //
    Bb = circulant_mat(b);


    // //construct matrix Hx
    //Aa = circulant_mat(a);
    //Bb = circulant_mat(b)
    // for(int m = 0; m < l; m++){
    //     vtemp = Aa.get_col(m);
    //     Hx.set_col(m, vtemp);
    // }
    // for(int o = l; o < 2*l; o++){
    //     vtemp = Bb.get_col(o-l);
    //     Hx.set_col(o, vtemp);
    // }
    // //find rank of Hx
    // rank = Hx.row_rank();


    //find roots of l
    //roots = NTL::FindRoots(l_pol);

    //find factors of l, a, and b
    l_factors = NTL::CanZass(l_pol, verbose);
    a_factors = NTL::CanZass(a_pol, verbose);
    b_factors = NTL::CanZass(b_pol, verbose);

    std::cout <<"l:"<< l_pol << "\n";
    std::cout <<"l_factors"<< l_factors << "\n";
    std::cout <<"a_pol:"<< a_pol << "\n";
    std::cout <<"a_factors"<< a_factors << "\n";
    std::cout <<"a:"<< a << "\n";
    std::cout <<"Aa:"<< Aa << "\n";
    std::cout <<"b_pol:"<< b_pol << "\n";
    std::cout <<"b_factors"<< b_factors << "\n";
    std::cout <<"b:"<< b << "\n";
    std::cout <<"Bb:"<< Bb << "\n";
    
    //find GCD of a and x^l-1
    g1 = NTL::GCD(a_pol, l_pol);
    g2 = NTL::GCD(g1, b_pol);

    //find degree of a
    g_deg = NTL::deg(g2);

    // d = common::quantum_dist_v2(Hx, Hx, 0);
    // Ad = common::classical_dist(Aa);
    
    m = NTL::MulMod(a_pol, b_pol, l_pol);

    c1_pol += NTL::GF2X(12,1) + NTL::GF2X(8,1)+ NTL::GF2X(4,1);
    c2_pol += 1;
    ab_pol += NTL::GF2X(14,1) + NTL::GF2X(13,1) + NTL::GF2X(4,1) +1;
    c_pol += NTL::GF2X(13,1) + NTL::GF2X(12,1) + NTL::GF2X(8,1) + NTL::GF2X(4,1);

    x1 = NTL::MulMod(b_pol, c1_pol, l_pol);
    x2 = NTL::MulMod(a_pol, c2_pol, l_pol);
    y = NTL::MulMod(ab_pol, c_pol, l_pol);
    //z = NTL::project(ab_pol, c_pol)
    
    // std::cout <<"Hx_rank:"<< rank << "\n";
    //std::cout <<"roots"<< roots << "\n";

    std::cout <<"x:"<< x1 <<", "<< x2 << "\n";
    std::cout <<"y:"<< y << "\n";
    //std::cout <<"z:"<< z << "\n";
    std::cout <<"GCD of a and l:"<< g1 << "\n";
    std::cout <<"GCD of a, b, and l:"<< g2 << "\n";
    std::cout <<"g_degree:"<< g_deg << "\n";

    //construct Hx
    for(int i = 0; i < l; i++){
        vtemp = Aa.get_col(i);
        Hx.set_col(i, vtemp);
    }
    for(int j = l; j < 2*l; j++){
        vtemp = Bb.get_col(j-l);
        Hx.set_col(j, vtemp);
    }

    std::cout << "Hx" << Hx << "\n";

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

    std::cout << "Hz" << Hz << "\n";

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

    //std::cout << "call qd" << "\n";
    //Hx = common::make_it_full_rank(Hx);
    //std::cout << "Hx" << Hx << "\n";
    //Hz = common::make_it_full_rank(Hz);
    d = common::quantum_dist_v2(Hz, Hx, 0);
    //itpp::GF2mat C = common::getC(Hx, Hz, 1);
    //itpp::bvec v2temp("1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0");
    //itpp::GF2mat C(1,26);
    //itpp::GF2mat C(1, 2*l);
    //C.set_row(0,v2temp);
    //std::cout << "C = " << C << "\n";
    //std::cout << "Hx*C = "<< Hx*C.transpose() << "\n";
    //itpp::GF2mat Q=Hz.concatenate_vertical(C);
    //std::cout << "rankHz: " << Hz.row_rank() << "," << Hz.rows() << "\n";
    //std::cout << "rankQ: " << Q.row_rank() <<"," <<Q.rows() << "\n";
    //std::cout << "Hx*Hz = " << Hx*Hz.transpose() << "\n";

    std::cout << "bicycle distance:" << d << "\n";
    // std::cout << "cyclic distance:" << Ad << "\n";


    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
    std::cout << "Run-time " << time_span.count() << " seconds.\n";
    std::cout << "Run-time " << time_span.count()/60 << " minutes.\n";
    std::cout << "Run-time " << time_span.count()/3600 << " hours.\n";

    return 0;
}
