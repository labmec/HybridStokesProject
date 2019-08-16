//
//  TPZSimulationData.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZSimulationData.h"

TPZSimulationData::TPZSimulationData()
{
    m_n_divs.resize(0);
    m_h_domain.resize(0);
    m_internal_order = 0;
    m_skeleton_order = 0;
    m_dimesion = 0;
    m_n_threads = 0;
    m_visco = 0;
    m_n_intrefs = 0.;
}

TPZSimulationData::~TPZSimulationData()
{
    
}

TPZSimulationData::TPZSimulationData(const TPZSimulationData & other)
{
    m_n_divs                           = other.m_n_divs;
    m_h_domain                         = other.m_h_domain;
    m_internal_order                   = other.m_internal_order;
    m_skeleton_order                   = other.m_skeleton_order;
    m_dimesion                         = other.m_dimesion;
    m_n_threads                        = other.m_n_threads;
    m_visco                            = other.m_visco;
    m_n_intrefs                        = other.m_n_intrefs;
}

TPZSimulationData & TPZSimulationData::operator=(const TPZSimulationData &other)
{
    if (this != & other) // prevent self-assignment
    {
        m_n_divs                           = other.m_n_divs;
        m_h_domain                         = other.m_h_domain;
        m_internal_order                   = other.m_internal_order;
        m_skeleton_order                   = other.m_skeleton_order;
        m_dimesion                         = other.m_dimesion;
        m_n_threads                        = other.m_n_threads;
        m_visco                            = other.m_visco;
        m_n_intrefs                        = other.m_n_intrefs;
    }
    return *this;
}

void TPZSimulationData::Print()
{
    
    std::cout << " TPZSimulationData class members : " << std::endl;
    std::cout << std::endl;
    std::cout << " m_n_divs = " << m_n_divs << std::endl;
    std::cout << " m_h_domain = " << m_h_domain << std::endl;
    std::cout << " m_internal_order = " << m_internal_order << std::endl;
    std::cout << " m_skeleton_order = " << m_skeleton_order << std::endl;
    std::cout << " m_dimesion = " << m_dimesion << std::endl;
    std::cout << " m_n_threads = " << m_n_threads << std::endl;
    std::cout << " m_visco = " << m_visco << std::endl;
    std::cout << " m_n_intrefs = " << m_n_intrefs << std::endl;
    std::cout << std::endl;
    
}

