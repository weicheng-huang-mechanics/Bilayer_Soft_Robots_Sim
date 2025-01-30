#ifndef SYMBOLICEQUATIONS_H
#define SYMBOLICEQUATIONS_H

#include <symengine/llvm_double.h>

using namespace SymEngine;


class symbolicEquations
{
public:
    symbolicEquations();

    void generatePenaltyFunctions();

    LLVMDoubleVisitor E_func;
    LLVMDoubleVisitor E_grad_dx_func;
    LLVMDoubleVisitor E_hess_dx_func;

    LLVMDoubleVisitor E_func_1;
    LLVMDoubleVisitor E_grad_dx_func_1;
    LLVMDoubleVisitor E_hess_dx_func_1;

    LLVMDoubleVisitor E_func_2;
    LLVMDoubleVisitor E_grad_dx_func_2;
    LLVMDoubleVisitor E_hess_dx_func_2;

    LLVMDoubleVisitor E_func_3;
    LLVMDoubleVisitor E_grad_dx_func_3;
    LLVMDoubleVisitor E_hess_dx_func_3;

    LLVMDoubleVisitor E_grad_dm_func;
    LLVMDoubleVisitor E_hess_dm_func;


    LLVMDoubleVisitor m1_1_func;
    LLVMDoubleVisitor m2_1_func;
    LLVMDoubleVisitor m1_2_func;
    LLVMDoubleVisitor m2_2_func;

    LLVMDoubleVisitor m1_1_func_1;
    LLVMDoubleVisitor m2_1_func_1;
    LLVMDoubleVisitor m1_2_func_1;
    LLVMDoubleVisitor m2_2_func_1;

    LLVMDoubleVisitor m1_1_func_2;
    LLVMDoubleVisitor m2_1_func_2;
    LLVMDoubleVisitor m1_2_func_2;
    LLVMDoubleVisitor m2_2_func_2;

    LLVMDoubleVisitor m1_1_func_3;
    LLVMDoubleVisitor m2_1_func_3;
    LLVMDoubleVisitor m1_2_func_3;
    LLVMDoubleVisitor m2_2_func_3;
    
    LLVMDoubleVisitor d1_1_func;

    LLVMDoubleVisitor E_m1_grad_func;
    LLVMDoubleVisitor E_m1_hess_func;

    LLVMDoubleVisitor E_m2_grad_func;
    LLVMDoubleVisitor E_m2_hess_func;

    LLVMDoubleVisitor E_m3_grad_func;
    LLVMDoubleVisitor E_m3_hess_func;

    LLVMDoubleVisitor E_trans1_grad_func;
    LLVMDoubleVisitor E_trans1_hess_func;
    
    LLVMDoubleVisitor E_trans2_grad_func;
    LLVMDoubleVisitor E_trans2_hess_func;

    LLVMDoubleVisitor E_trans3_grad_func;
    LLVMDoubleVisitor E_trans3_hess_func;

    LLVMDoubleVisitor E_parallel_grad_func;
    LLVMDoubleVisitor E_parallel_hess_func;

private:
    bool symbolic_cse;
    int opt_level;

    // Helper functions for symbolic differentiation process
    void subtract_matrix(const DenseMatrix &A, const DenseMatrix &B, DenseMatrix &C);
    void get_norm(const DenseMatrix &num, RCP<const Basic> &C);
    void convert_to_unit_vector(const DenseMatrix &num, DenseMatrix &C);
    void cross_product(const DenseMatrix &u, const DenseMatrix &v, DenseMatrix &c);
    void parallelTransport(const DenseMatrix &d1_1, const DenseMatrix &t1, const DenseMatrix &t2, DenseMatrix &d1_2);

    RCP<const Basic> x1s_x;
    RCP<const Basic> x1s_y;
    RCP<const Basic> x1s_z;
    RCP<const Basic> x1e_x;
    RCP<const Basic> x1e_y;
    RCP<const Basic> x1e_z;
    RCP<const Basic> x2s_x;
    RCP<const Basic> x2s_y;
    RCP<const Basic> x2s_z;
    RCP<const Basic> x2e_x;
    RCP<const Basic> x2e_y;
    RCP<const Basic> x2e_z;
    RCP<const Basic> K;
    RCP<const Basic> theta_1;
    RCP<const Basic> theta_2;

    RCP<const Basic> x1s_x0;
    RCP<const Basic> x1s_y0;
    RCP<const Basic> x1s_z0;
    RCP<const Basic> x1e_x0;
    RCP<const Basic> x1e_y0;
    RCP<const Basic> x1e_z0;
    RCP<const Basic> x2s_x0;
    RCP<const Basic> x2s_y0;
    RCP<const Basic> x2s_z0;
    RCP<const Basic> x2e_x0;
    RCP<const Basic> x2e_y0;
    RCP<const Basic> x2e_z0;


    RCP<const Basic> h1;

    RCP<const Basic> d1_1x0;
    RCP<const Basic> d1_1y0;
    RCP<const Basic> d1_1z0;
    RCP<const Basic> d1_2x0;
    RCP<const Basic> d1_2y0;
    RCP<const Basic> d1_2z0;

    RCP<const Basic> d2_1x0;
    RCP<const Basic> d2_1y0;
    RCP<const Basic> d2_1z0;
    RCP<const Basic> d2_2x0;
    RCP<const Basic> d2_2y0;
    RCP<const Basic> d2_2z0;
    
    RCP<const Basic> d1_1x;
    RCP<const Basic> d1_1y;
    RCP<const Basic> d1_1z;
    RCP<const Basic> d1_2x;
    RCP<const Basic> d1_2y;
    RCP<const Basic> d1_2z;

    RCP<const Basic> d2_1x;
    RCP<const Basic> d2_1y;
    RCP<const Basic> d2_1z;
    RCP<const Basic> d2_2x;
    RCP<const Basic> d2_2y;
    RCP<const Basic> d2_2z;


    RCP<const Basic> m1_x;
    RCP<const Basic> m1_y;
    RCP<const Basic> m1_z;

    RCP<const Basic> m2_x;
    RCP<const Basic> m2_y;
    RCP<const Basic> m2_z;
};

#endif