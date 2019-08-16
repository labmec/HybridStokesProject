//
//  TPZSimulationData.h
//  PZ
//
//  Created by Omar on 8/28/18.
//
//

#ifndef TPZSimulationData_h
#define TPZSimulationData_h

#include <stdio.h>
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"


/** @brief Object conatining several kind of informations being used anytime and anywhere */
class TPZSimulationData
{
    
protected:
    
    /** Number of coarse elements in each axis*/
    TPZVec<int> m_n_divs;
    
    /** Number of refinements in each internal element*/
    int m_n_intrefs;
    
    /** Size of domain*/
    TPZVec<REAL> m_h_domain;
    
    /** Polynomial order for internal elements */
    int m_internal_order;
    
    /** Polynomial order for internal elements */
    int m_skeleton_order;
    
    /** Physical dimension of the domain */
    int m_dimesion;
    
    /** Number of thread */
    int m_n_threads;

    /** Viscosity coeficient */
    REAL m_visco;
    
public:
    
    /** default constructor */
    TPZSimulationData();
    
    /** default constructor */
    TPZSimulationData(const TPZSimulationData & other);
    
    /** default constructor */
    TPZSimulationData &operator=(const TPZSimulationData &other);
    
    /** destructor */
    ~TPZSimulationData();
    
    /** Print object attributes */
    void Print();

    /** Set the number of coarse elements in each axis*/
    void SetCoarseDivisions(TPZVec<int> n_divs){
        m_n_divs = n_divs;
    }

    /** Get normal stiffness for each fracture*/
    TPZVec<int> GetCoarseDivisions(){
        return m_n_divs;
    }

    /** Set the number of refinements in each internal element*/
    void SetNInterRefs(int n_refs){
        m_n_intrefs = n_refs;
    }
    
    /** Get the number of refinements in each internal element*/
    int GetNInterRefs(){
        return m_n_intrefs;
    }
    
    /** Set the size of domain*/
    void SetDomainSize(TPZVec<REAL> h_domain){
        m_h_domain = h_domain;
    }

    /** Get the size of domain*/
    TPZVec<REAL> GetDomainSize(){
        return m_h_domain;
    }
    
    /** Set polynomial order for internal elements */
    void SetInternalOrder(int internal_order){
        m_internal_order = internal_order;
    }

    /** Get polynomial order for internal elements */
    int GetInternalOrder(){
        return m_internal_order;
    }
    
    /** Set polynomial order for skeleton elements */
    void SetSkeletonOrder(int skeleton_order){
        m_skeleton_order = skeleton_order;
    }
    
    /** Get polynomial order for skeleton elements */
    int GetSkeletonOrder(){
        return m_skeleton_order;
    }

    /** Set the physical dimension of the domain */
    void SetDimension(int dim){
        m_dimesion = dim;
    }
    
    /** Get the physical dimension of the domain */
    int GetDimension(){
        return m_dimesion;
    }

    /** Set the number of threads */
    void SetNthreads(int nthreads){
        m_n_threads = nthreads;
    }
    
    /** Get the number of threads */
    int GetNthreads(){
        return m_n_threads;
    }

    /** Set the viscosity coeficient */
    void SetViscosity(REAL viscosity){
        m_visco = viscosity;
    }
    
    /** Get the viscosity coeficient */
    REAL GetViscosity(){
        return m_visco;
    }
    
};

#endif /* TPZSimulationData_h */
