//
//  manipContacts.h
//  
//
//  Created by David Hodgson on 06/07/2021.
//

#ifndef manipContacts_h
#define manipContacts_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include "ascent/Ascent.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]
using namespace Rcpp;
using namespace std;
using namespace Eigen;
using namespace asc;

class RunInterventions2 {
    
    int A;

    vector<double > p_mat;
    vector<vector<double > > cnt_matrix_p;
    vector<vector<double > > cnt_matrix_p_h;
    vector<vector<double > > pwp_p;
    vector<vector<double > > pwn_p;
    vector<vector<double > > nwp_p;
    vector<vector<double > > nwn_p;
    
    vector<vector<double > > cnt_matrix_c;
    vector<vector<double > > cnt_matrix_c_h;
    vector<vector<double > > pwp_c;
    vector<vector<double > > pwn_c;
    vector<vector<double > > nwp_c;
    vector<vector<double > > nwn_c;
        
    /**********************************************************************/
    /** CwX Contacts for persons who are cocooned ****/
    /**********************************************************************/
    
    vector<vector<double > > get_cwn(double prop_c, char s)
    {
        vector<vector<double > > cwn_e(A,vector<double>(A, 0));
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (nwn_p[i][j]);
            
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwn_e[i][j] = (nwn_p[i][j])*(1-prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (pwn_p[i][j]);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwn_e[i][j] = (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(1-prop_c);

        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (nwn_c[i][j]);
            
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwn_e[i][j] = (nwn_c[i][j])*(1-prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (pwn_c[i][j]);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwn_e[i][j] = (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return cwn_e;
    }
    
    vector<vector<double > > get_cwp(double prop_c, char s)
    {
        vector<vector<double > > cwp_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwp_e[i][j] = (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(p_mat[j])*(1-prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwp_e[i][j] = (pwp_p[i][j])*(1-prop_c); // parental w. infant plus non-parental with infant who are cocooned
    
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwp_e[i][j] = (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(p_mat[j])*(1-prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwp_e[i][j] = (pwp_c[i][j])*(1-prop_c); // parental w. infant plus non-parental with infant who are cocooned
        }
        else
            cout << "Error cont" << endl;
        
        return cwp_e;
    }
    
    vector<vector<double > > get_cwc(double prop_c, char s)
    {
        vector<vector<double > > cwc_e(A,vector<double>(A,0));

        if (s == 'p')
        {
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwc_e[i][j] = (nwn_p[i][j])*(prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwc_e[i][j] = cnt_matrix_p_h[i][j] + (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(prop_c);
            
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwc_e[i][j] = cnt_matrix_p_h[i][j]/2.0 + (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(p_mat[j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwc_e[i][j] = (pwp_p[i][j])*(prop_c); // parental w. infant plus non-parental with infant who are cocooned
    
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwc_e[i][j] = (nwn_c[i][j])*(prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwc_e[i][j] = cnt_matrix_c_h[i][j] + (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(prop_c);
            
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwc_e[i][j] = cnt_matrix_c_h[i][j]/2.0 + (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(p_mat[j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwc_e[i][j] = (pwp_c[i][j])*(prop_c); // parental w. infant plus non-parental with infant who are cocooned

        }
        else
            cout << "Error cont" << endl;
        return cwc_e;
    }
    
    /**********************************************************************/
    /** PwX Contacts for persons who are mothers but not cocooned ****/
    /**********************************************************************/
    
    
    vector<vector<double > > get_pwn(double prop_c, char s)
    {
        vector<vector<double > > pwn_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < A; j++) // all
                    pwn_e[i][j] = pwn_p[i][j]; // all contacts
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwn_e[i][j] = cnt_matrix_p_h[i][j] + (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(1-prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < A; j++) // all
                    pwn_e[i][j] = pwn_c[i][j]; // all contacts
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwn_e[i][j] = cnt_matrix_c_h[i][j] + (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return pwn_e;
    }
    
    vector<vector<double > > get_pwp(double prop_c, char s)
    {
        vector<vector<double > > pwp_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // all
                    pwp_e[i][j] = pwp_p[i][j]*(1-prop_c); // all contacts
        }
        else if (s == 'c')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // all
                    pwp_e[i][j] = pwp_c[i][j]*(1-prop_c); // all contacts
        }
        else
            cout << "Error cont" << endl;
        
        return pwp_e;
    }
    
    vector<vector<double > > get_pwc(double prop_c, char s)
    {
        // Non cocooned parent -> cocooned infant (family and outside family)
        vector<vector<double > > pwc_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwc_e[i][j] = (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    pwc_e[i][j] = pwn_p[i][j]*prop_c; //
        }
        else if (s == 'c')
        {
            // Parent -> parents (all)
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwc_e[i][j] = (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    pwc_e[i][j] = pwn_c[i][j]*prop_c; //
        }
        else
            cout << "Error cont" << endl;
        
        return pwc_e;
    }
    
    
 /**********************************************************************/
/** NwX Contacts for persons who are neither cocooned nor mothers ****/
 /**********************************************************************/
    
    vector<vector<double > > get_nwn(double prop_c, char s)
    {
        // cnt_matrix_p cnt_matrix_p_h
        vector<vector<double > > nwn_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < A; j++) //     infants
                    nwn_e[i][j] = nwn_p[i][j]; //number of contacts made by population with non-cocooned prop_ci
            
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwn_e[i][j] = nwn_p[i][j]*(1-prop_c); //number
        }
        else if (s == 'c')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < A; j++) //     infants
                    nwn_e[i][j] = nwn_c[i][j]; //number of contacts made by population with non-cocooned prop_ci
            
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwn_e[i][j] = nwn_c[i][j]*(1-prop_c); //numbe

        }
        else
            cout << "Error cont" << endl;

        return nwn_e;
    }
    
    vector<vector<double > > get_nwp(double prop_c, char s)
    {
        // Population with non-cocooned parents)
        vector<vector<double > > nwp_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = (cnt_matrix_p_h[i][j]/2.0) + (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*p_mat[j]*(1-prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = cnt_matrix_p[i][j]*p_mat[j]*(1-prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = (cnt_matrix_c_h[i][j]/2.0) + (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*p_mat[j]*(1-prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = cnt_matrix_c[i][j]*p_mat[j]*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return nwp_e;
    }
    
    vector<vector<double > > get_nwc(double prop_c, char s)
    {
        
        vector<vector<double > > nwc_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwc_e[i][j] = nwn_p[i][j]*(prop_c); //numbe
            
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*p_mat[j]*(prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = cnt_matrix_p[i][j]*p_mat[j]*(prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwc_e[i][j] = nwn_c[i][j]*(prop_c); //numbe
            
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*p_mat[j]*(prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = cnt_matrix_c[i][j]*p_mat[j]*(prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return nwc_e;
    }
};

#endif /* manipContacts_h */
