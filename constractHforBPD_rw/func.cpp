#include<iostream>
#include<cmath>
#include<itpp/itbase.h>
#include<NTL/GF2X.h>
#include<fstream>
#include "weilei_lib/weilei_lib.h"


// write a number into binary form
itpp::bvec binary_num(int r1, int x1){

    itpp::bvec co(r1);
    int j1 = r1 - 1;

    co.zeros();

    while (x1 > 0) { 
            co(j1) = x1 % 2; 
            x1 = x1 / 2; 
            j1--; 
    }
    
    return co;

}

int weight(itpp::bvec a2){
    int n2 = a2.size();
    int w2 = 0;
    for (int i2 = 0; i2 < n2; i2++) {
        if (a2(i2) == 1){
            w2++;
        }
    }

    return w2;
}

itpp::GF2mat circulant_mat(itpp::bvec a3){
    int l3 = a3.size();
    itpp::bvec a3_temp(l3),a3_temp2(l3);
    itpp::GF2mat Aa3(l3,l3);
    for(int j3 = 0; j3 < l3; j3++){
        a3_temp(j3) = a3(l3-j3-1);    
    }
    Aa3.set_row(l3-1,a3_temp);

    for(int k3 = l3-2; k3 >= 0; k3--){
        for(int m3 = 0; m3 < l3; m3++){
            a3_temp2(m3) = a3_temp((m3+1)%l3);
        }
        a3_temp = a3_temp2;
        Aa3.set_row(k3, a3_temp);
    }

    return Aa3;
}

void set_weight(itpp::bvec a4, int begin4, int w4, int l4){

    int cof, k4, d4, rate4;
    long coff;
    NTL::GF2X a4_pol, g4;
    NTL::GF2X l4_pol = NTL::GF2X(l4,1) + 1;
    itpp::GF2mat Aa4(l4,l4), TAa4(l4,l4), Hx4(l4,2*l4);
    itpp::bvec vtemp(l4);

    //std::ofstream data;
    //data.open ("data.txt");
    
    if(w4 > 1){
        for(int i4 = begin4; i4 <= l4-w4+1; i4++){
            a4(i4) = 1;
            set_weight(a4, i4+1, w4-1,l4);
            a4(i4) = 0;
        }
    }
    else{

        for(int j4 = 0; j4 < l4; j4++){
            cof = a4(j4);
            coff = cof;
            a4_pol += NTL::GF2X(j4,coff);
        }
        
        g4 = GCD(a4_pol, l4_pol);

        if(g4 != 1){
        
            k4 = 2 * NTL::deg(g4);
            if(k4 != 0){

                Aa4 = circulant_mat(a4);
                TAa4 = Aa4.transpose();

                for(int m4 = 0; m4 < l4; m4++){
                    vtemp = Aa4.get_col(m4);
                    Hx4.set_col(m4, vtemp);
                }

                for(int o4 = l4; o4 < 2*l4; o4++){
                    vtemp = TAa4.get_col(o4-l4);
                    Hx4.set_col(o4, vtemp);
                }
                
                d4 = common::quantum_dist_v2(Hx4, Hx4, 0);

                if(d4 > 2 & d4 < 999){
                    //rate4 = k4/l4/2;
                    //count4 = count4 + 1;
                    std::cout << "[" << 2*l4 <<", "<< k4 <<", "
                        << d4 << "]" << "non zero digit of a: ";
                    for(int p4 = 0; p4 < l4; p4++){
                        if(a4(p4)!= 0){
                            std::cout << p4 << ", ";
                        }
                    }
                    std::cout << "\n"; 
                }
                
                
            }
        }

        a4_pol = 0;
    }


}
    

