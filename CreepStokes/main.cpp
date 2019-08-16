

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "StokesTest.h"
#include "HybridBrinkmanTest.h"
#include "TPZStokesMaterial.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include "pztrnsform.h"

#define TEST_DOMAINS
//#define APP_CURVE

const int SpaceHDiv = 1; //Velocidade em subespaço de H(div)
const int SpaceContinuous = 2; //Velocidade em subespaço de [H1]ˆ2
const int SpaceDiscontinuous = 3; //Velociadade em subespaço de H(Ph) - Ph: partição
const REAL Pi=M_PI;

//Verificação dos modelos:
#ifdef TEST_DOMAINS

const REAL visco=1., permeability=1., theta=-1.; //Coeficientes: viscosidade, permeabilidade, fator simetria

bool StokesDomain = false;
bool HybridBrinkmanDomain = true;

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e12;
//    gRefDBase.InitializeAllUniformRefPatterns();
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    
    REAL hx=2.,hy=2.; //Dimensões em x e y do domínio
    //double hx=Pi,hy=2.;
    int h_level = 0;
    int nx=h_level+1 ,ny=h_level+1; //Número de nos em x  y
    int pOrder = 0; //Ordem polinomial de aproximação
    
    TPZVec<REAL> h_s(3,0);
    h_s[0]=2.,h_s[1]=2.,h_s[2]=2.; //Dimensões em x e y do domínio

    if (HybridBrinkmanDomain){
        
        int pOrder = 2;
        
        HDivPiola = 1;
        for (int it=0; it<=0; it++) {
            //h_level = pow(2., 2+it);
            h_level = 2;
            
            TPZVec<int> n_s(3,0.);
            n_s[0]=h_level ,n_s[1]=h_level;
            
            n_s[2]=h_level; //Obs!!
            
            REAL visc = 1.0; //->Darcy
            
            HybridBrinkmanTest  * Test2 = new HybridBrinkmanTest();
            //Test2->Set3Dmesh();
            //Test2->SetElType(ETetraedro);
            //Test2->SetHdivPlus();

            TPZTransform<STATE> Transf(3,3), InvTransf(3,3);
            Test2->SetTransform(Transf, InvTransf);

            REAL rot_x = 5.;
            REAL rot_z = 44.;
            REAL rot_y = -85.;
            rot_z = rot_z*Pi/180.;
            rot_y = rot_y*Pi/180.;
            rot_z = rot_z*Pi/180.;
            
            //Test2->SetRotation3DMatrix(rot_x,rot_y,rot_z);
            //Test2->SetAllRefine();
            Test2->Run(SpaceHDiv, pOrder, n_s, h_s,visc);
            
        }
        
    }else if(StokesDomain)
    {
        pOrder = 3;

        TPZVec<STATE> S0(13,0.);
        S0[0]=0.0000001,S0[1]=1.,S0[2]=3.,S0[3]=5.,S0[4]=10.,S0[5]=15.,S0[6]=20.,S0[7]=25.,S0[8]=30.,S0[9]=35.,S0[10]=40.,S0[11]=45.,S0[12]=50.;
        HDivPiola = 0;
        
        hx=2., hy=2.;

        for (int it=0; it<=12; it++) {
            //h_level = pow(2., 2+it);
            h_level = 8;
            //Coeficiente estabilização (Stokes)
            STATE hE=hx/h_level;
            STATE s0=S0[it];
            STATE sigma=s0*(pOrder*pOrder)/hE;


            nx=h_level+1 ,ny=h_level+1;
            hE=hx/h_level;
            sigma=s0*(pOrder*pOrder)/hE;
            StokesTest  * Test1 = new StokesTest();
            Test1->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visco,theta,sigma);
            //h_level = h_level*2;
        }
        
    }
    
    return 0;
}
#endif




