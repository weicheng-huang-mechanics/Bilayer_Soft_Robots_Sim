#include "symbolicEquations.h"
#include <chrono>


symbolicEquations::symbolicEquations() {
    // nodes
    x1s_x = symbol("x1s_x");
    x1s_y = symbol("x1s_y");
    x1s_z = symbol("x1s_z");
    x1e_x = symbol("x1e_x");
    x1e_y = symbol("x1e_y");
    x1e_z = symbol("x1e_z");
    theta_1 = symbol("theta_1");

    x2s_x = symbol("x2s_x");
    x2s_y = symbol("x2s_y");
    x2s_z = symbol("x2s_z");
    x2e_x = symbol("x2e_x");
    x2e_y = symbol("x2e_y");
    x2e_z = symbol("x2e_z");
    theta_2 = symbol("theta_2");

    x1s_x0 = symbol("x1s_x0");
    x1s_y0 = symbol("x1s_y0");
    x1s_z0 = symbol("x1s_z0");
    x1e_x0 = symbol("x1e_x0");
    x1e_y0 = symbol("x1e_y0");
    x1e_z0 = symbol("x1e_z0");

    x2s_x0 = symbol("x2s_x0");
    x2s_y0 = symbol("x2s_y0");
    x2s_z0 = symbol("x2s_z0");
    x2e_x0 = symbol("x2e_x0");
    x2e_y0 = symbol("x2e_y0");
    x2e_z0 = symbol("x2e_z0");



    // stiffness    
    K = symbol("K");
    h1 = symbol("h1");

    d1_1x0 = symbol("d1_1x0");
    d1_1y0 = symbol("d1_1y0");
    d1_1z0 = symbol("d1_1z0");

    d1_2x0 = symbol("d1_2x0");
    d1_2y0 = symbol("d1_2y0");
    d1_2z0 = symbol("d1_2z0");

    d2_1x0 = symbol("d2_1x0");
    d2_1y0 = symbol("d2_1y0");
    d2_1z0 = symbol("d2_1z0");

    d2_2x0 = symbol("d2_2x0");
    d2_2y0 = symbol("d2_2y0");
    d2_2z0 = symbol("d2_2z0");

    d1_1x = symbol("d1_1x");
    d1_1y = symbol("d1_1y");
    d1_1z = symbol("d1_1z");

    d1_2x = symbol("d1_2x");
    d1_2y = symbol("d1_2y");
    d1_2z = symbol("d1_2z");

    d2_1x = symbol("d2_1x");
    d2_1y = symbol("d2_1y");
    d2_1z = symbol("d2_1z");

    d2_2x = symbol("d2_2x");
    d2_2y = symbol("d2_2y");
    d2_2z = symbol("d2_2z");

    symbolic_cse = true;
    opt_level = 3;
}


// For some reason SymEngine doesn't have this implemented X_X
void symbolicEquations::subtract_matrix(const DenseMatrix &A, const DenseMatrix &B, DenseMatrix &C) {
    assert((A.nrows() == B.nrows()) && (A.ncols() == B.ncols()));
    for (unsigned i=0; i < A.nrows(); i++) {
        for (unsigned j=0; j < A.ncols(); j++) {
            C.set(i, j, sub(A.get(i, j), B.get(i, j)));
        }
    }
}

void symbolicEquations::cross_product(const DenseMatrix &u, const DenseMatrix &v, DenseMatrix &c) {
    // Ensure u and v are 3x1 matrices
    if (u.nrows() != 3 || u.ncols() != 1 || v.nrows() != 3 || v.ncols() != 1) {
        throw std::invalid_argument("Both input matrices must be 3x1 vectors.");
    }

    // Resize the output matrix to 3x1
    c.resize(3, 1);

    // Compute cross product components
    c.set(0, 0, sub(mul(u.get(1, 0), v.get(2, 0)), mul(u.get(2, 0), v.get(1, 0)))); // cx = u2*v3 - u3*v2
    c.set(1, 0, sub(mul(u.get(2, 0), v.get(0, 0)), mul(u.get(0, 0), v.get(2, 0)))); // cy = u3*v1 - u1*v3
    c.set(2, 0, sub(mul(u.get(0, 0), v.get(1, 0)), mul(u.get(1, 0), v.get(0, 0)))); // cz = u1*v2 - u2*v1
}



void symbolicEquations::get_norm(const DenseMatrix &num, RCP<const Basic> &C) {
    DenseMatrix tmp(num.nrows(), num.ncols());
    num.elementwise_mul_matrix(num, tmp);
    C = sqrt(add(tmp.as_vec_basic()));
}


void symbolicEquations::convert_to_unit_vector(const DenseMatrix &num, DenseMatrix &C) {
    DenseMatrix tmp(num.nrows(), num.ncols());
    num.elementwise_mul_matrix(num, tmp);
    auto norm = sqrt(add(tmp.as_vec_basic()));
    for (unsigned i=0; i < num.nrows(); i++) {
        for (unsigned j=0; j < num.ncols(); j++) {
            C.set(i, j, div(num.get(i, j), norm));
        }
    }
}

void symbolicEquations::parallelTransport(const DenseMatrix &d1_1, const DenseMatrix &t1, const DenseMatrix &t2, DenseMatrix &d1_2) 
{
    DenseMatrix b(3, 1);
    DenseMatrix n1(3, 1);
    DenseMatrix n2(3, 1);

    cross(t1, t2, b); // b = t1 x t2

    convert_to_unit_vector(b, b); // b = b/norm(b)

    cross(t1, b, n1); // n1 = t1 x b
    cross(t2, b, n2); // n2 = t2 x b
 
    // d1_2=d1_1.dot(t1)*t2+d1_1.dot(n1)*n2+d1_1.dot(b)*b;

    DenseMatrix temp1(3, 1);
    DenseMatrix temp2(3, 1);
    DenseMatrix temp3(3, 1);

    d1_1.elementwise_mul_matrix(t1, temp1); //temp1 = d1_1 * t1
    d1_1.elementwise_mul_matrix(n1, temp2); //temp2 = d1_1 * n1
    d1_1.elementwise_mul_matrix(b, temp3); //temp3 = d1_1 * b

    t2.mul_scalar(add(temp1.as_vec_basic()), temp1);
    n2.mul_scalar(add(temp2.as_vec_basic()), temp2);
    b.mul_scalar(add(temp3.as_vec_basic()), temp3);

    temp1.add_matrix(temp2, d1_2);
    d1_2.add_matrix(temp3, d1_2);

    d1_2.elementwise_mul_matrix(t2, temp1);
    t2.mul_scalar(add(temp1.as_vec_basic()), temp1);
    subtract_matrix(d1_2, temp1, d1_2);
    convert_to_unit_vector(d1_2, d1_2);    
}

void symbolicEquations::generatePenaltyFunctions(){
    DenseMatrix x1s({x1s_x, x1s_y, x1s_z});
    DenseMatrix x1e({x1e_x, x1e_y, x1e_z});
    DenseMatrix x2s({x2s_x, x2s_y, x2s_z});
    DenseMatrix x2e({x2e_x, x2e_y, x2e_z});

    DenseMatrix x1s0({x1s_x0, x1s_y0, x1s_z0});
    DenseMatrix x1e0({x1e_x0, x1e_y0, x1e_z0});
    DenseMatrix x2s0({x2s_x0, x2s_y0, x2s_z0});
    DenseMatrix x2e0({x2e_x0, x2e_y0, x2e_z0});

    DenseMatrix d1_10({d1_1x0, d1_1y0, d1_1z0});
    DenseMatrix d1_20({d1_2x0, d1_2y0, d1_2z0});
    DenseMatrix d2_10({d2_1x0, d2_1y0, d2_1z0});
    DenseMatrix d2_20({d2_2x0, d2_2y0, d2_2z0});

    vec_basic func_inputs{x1s_x, x1s_y, x1s_z, theta_1,
                        x1e_x, x1e_y, x1e_z,
                        x2s_x, x2s_y, x2s_z, theta_2,
                        x2e_x, x2e_y, x2e_z,
                        x1s_x0, x1s_y0, x1s_z0,
                        x1e_x0, x1e_y0, x1e_z0,
                        x2s_x0, x2s_y0, x2s_z0,
                        x2e_x0, x2e_y0, x2e_z0,
                        d1_1x0, d1_1y0, d1_1z0,
                        d1_2x0, d1_2y0, d1_2z0,
                        d2_1x0, d2_1y0, d2_1z0,
                        d2_2x0, d2_2y0, d2_2z0};
                         
    func_inputs.push_back(K);
    func_inputs.push_back(h1);

    DenseMatrix m1(3, 1);
    DenseMatrix m2(3, 1);

    DenseMatrix d1_1(3, 1);
    DenseMatrix d1_2(3, 1);

    DenseMatrix d2_1(3, 1);
    DenseMatrix d2_2(3, 1);

    // compute the tangent
    DenseMatrix t1(3, 1);
    DenseMatrix t10(3, 1);

    DenseMatrix t2(3, 1);
    DenseMatrix t20(3, 1);

    DenseMatrix m1_1(3, 1);
    DenseMatrix m1_2(3, 1);

    DenseMatrix m2_1(3, 1);
    DenseMatrix m2_2(3, 1);

    subtract_matrix(x1e, x1s, t1);
    convert_to_unit_vector(t1, t1);

    subtract_matrix(x1e0, x1s0, t10);
    convert_to_unit_vector(t10, t10);

    subtract_matrix(x2e, x2s, t2);
    convert_to_unit_vector(t2, t2);

    subtract_matrix(x2e0, x2s0, t20);
    convert_to_unit_vector(t20, t20);

    parallelTransport(d1_10, t10, t1, d1_1);
    cross(t1, d1_1, d1_2);

    parallelTransport(d2_10, t20, t2, d2_1);
    cross(t2, d2_1, d2_2);
    d1_1_func.init(func_inputs, d2_1.as_vec_basic(), symbolic_cse, opt_level);


    DenseMatrix temp1(3, 1);
    DenseMatrix temp2(3, 1);

    d1_1.mul_scalar(cos(theta_1), temp1);
    d1_2.mul_scalar(sin(theta_1), temp2);
    temp1.add_matrix(temp2, m1_1);

    d1_1.mul_scalar(sin(theta_1), temp1);
    d1_2.mul_scalar(cos(theta_1), temp2);
    subtract_matrix(temp2, temp1, m1_2);

    d2_1.mul_scalar(cos(theta_2), temp1);
    d2_2.mul_scalar(sin(theta_2), temp2);
    temp1.add_matrix(temp2, m2_1);

    d2_1.mul_scalar(sin(theta_2), temp1);
    d2_2.mul_scalar(cos(theta_2), temp2);
    subtract_matrix(temp2, temp1, m2_2);

    m1_1_func.init(func_inputs, m1_1.as_vec_basic(), symbolic_cse, opt_level);
    m2_1_func.init(func_inputs, m2_1.as_vec_basic(), symbolic_cse, opt_level);
    m1_2_func.init(func_inputs, m1_2.as_vec_basic(), symbolic_cse, opt_level);
    m2_2_func.init(func_inputs, m2_2.as_vec_basic(), symbolic_cse, opt_level);
    
    // define parallel item
    m1_1.elementwise_mul_matrix(m2_1, temp1);
    RCP <const Basic> parallel1 = pow(sub(add(temp1.as_vec_basic()), one), 2);    
    m1_2.elementwise_mul_matrix(m2_2, temp1);
    RCP<const Basic> parallel2 = pow(sub(add(temp1.as_vec_basic()), one), 2);
    t1.elementwise_mul_matrix(t2, temp1);
    RCP<const Basic> parallel3 = pow(sub(add(temp1.as_vec_basic()), one), 2);
    RCP<const Basic> parallel = add(add(parallel1, parallel2), parallel3);
    
    // define the penalty energy
    DenseMatrix x1(3, 1);
    x1s.add_matrix(x1e, x1);
    x1.mul_scalar(div(one, two), x1);

    DenseMatrix x2(3, 1);
    x2s.add_matrix(x2e, x2);
    x2.mul_scalar(div(one, two), x2);

    DenseMatrix dp(3, 1);
    subtract_matrix(x1, x2, dp);

    dp.elementwise_mul_matrix(m1_1, temp1);
    RCP<const Basic> trans1 = pow(sub(add(temp1.as_vec_basic()), h1), 2);
    dp.elementwise_mul_matrix(m1_2, temp1);
    RCP<const Basic> trans2 = pow(add(temp1.as_vec_basic()), 2);
    dp.elementwise_mul_matrix(t1, temp1);
    RCP<const Basic> trans3 = pow(add(temp1.as_vec_basic()), 2);
    RCP<const Basic> trans = add(add(trans1, trans2), trans3);

    // RCP<const Basic> E = mul(div(one, two), mul(K, add(trans, parallel)));  
    RCP<const Basic> E = mul(div(one, two), mul(K, trans));
    vec_basic nodes_vec {x1s_x, x1s_y, x1s_z, theta_1,
                         x1e_x, x1e_y, x1e_z,
                         x2s_x, x2s_y, x2s_z, theta_2,
                         x2e_x, x2e_y, x2e_z};
    DenseMatrix nodes {nodes_vec};


    DenseMatrix E_potential{{E}};
    DenseMatrix E_gradient(1, 14);
    DenseMatrix E_hessian(14, 14);
    jacobian(E_potential, nodes, E_gradient);    
    jacobian(E_gradient, nodes, E_hessian);


    auto start = std::chrono::high_resolution_clock::now();
    E_func.init(func_inputs, E_potential.as_vec_basic(), symbolic_cse, opt_level);
    E_grad_dx_func.init(func_inputs, E_gradient.as_vec_basic(), symbolic_cse, opt_level);
    E_hess_dx_func.init(func_inputs, E_hessian.as_vec_basic(), symbolic_cse, opt_level);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);



    // write func for d1 cross product is zero
    d1_1 = d1_10;
    cross(t1, d1_1, d1_2);

    parallelTransport(d2_10, t20, t2, d2_1);
    cross(t2, d2_1, d2_2);

    d1_1.mul_scalar(cos(theta_1), temp1);
    d1_2.mul_scalar(sin(theta_1), temp2);
    temp1.add_matrix(temp2, m1_1);

    d1_1.mul_scalar(sin(theta_1), temp1);
    d1_2.mul_scalar(cos(theta_1), temp2);
    subtract_matrix(temp2, temp1, m1_2);

    d2_1.mul_scalar(cos(theta_2), temp1);
    d2_2.mul_scalar(sin(theta_2), temp2);
    temp1.add_matrix(temp2, m2_1);

    d2_1.mul_scalar(sin(theta_2), temp1);
    d2_2.mul_scalar(cos(theta_2), temp2);
    subtract_matrix(temp2, temp1, m2_2);

    m1_1_func_1.init(func_inputs, m1_1.as_vec_basic(), symbolic_cse, opt_level);
    m2_1_func_1.init(func_inputs, m2_1.as_vec_basic(), symbolic_cse, opt_level);
    m1_2_func_1.init(func_inputs, m1_2.as_vec_basic(), symbolic_cse, opt_level);
    m2_2_func_1.init(func_inputs, m2_2.as_vec_basic(), symbolic_cse, opt_level);
    
    // define parallel item
    m1_1.elementwise_mul_matrix(m2_1, temp1);
    parallel1 = pow(sub(add(temp1.as_vec_basic()), one), 2);    
    m1_2.elementwise_mul_matrix(m2_2, temp1);
    parallel2 = pow(sub(add(temp1.as_vec_basic()), one), 2);
    t1.elementwise_mul_matrix(t2, temp1);
    parallel3 = pow(sub(add(temp1.as_vec_basic()), one), 2);
    parallel = add(add(parallel1, parallel2), parallel3);
    
    // define the penalty energy
    x1s.add_matrix(x1e, x1);
    x1.mul_scalar(div(one, two), x1);

    x2s.add_matrix(x2e, x2);
    x2.mul_scalar(div(one, two), x2);

    subtract_matrix(x1, x2, dp);

    dp.elementwise_mul_matrix(m1_1, temp1);
    trans1 = pow(sub(add(temp1.as_vec_basic()), h1), 2);
    dp.elementwise_mul_matrix(m1_2, temp1);
    trans2 = pow(add(temp1.as_vec_basic()), 2);
    dp.elementwise_mul_matrix(t1, temp1);
    trans3 = pow(add(temp1.as_vec_basic()), 2);
    trans = add(add(trans1, trans2), trans3);

    // E = mul(div(one, two), mul(K, add(trans, parallel)));
    E = mul(div(one, two), mul(K, trans));

    DenseMatrix E_potential1{{E}};
    DenseMatrix E_gradient1(1, 14);
    DenseMatrix E_hessian1(14, 14);
    jacobian(E_potential1, nodes, E_gradient1);    
    jacobian(E_gradient1, nodes, E_hessian1);

    E_func_1.init(func_inputs, E_potential1.as_vec_basic(), symbolic_cse, opt_level);
    E_grad_dx_func_1.init(func_inputs, E_gradient1.as_vec_basic(), symbolic_cse, opt_level);
    E_hess_dx_func_1.init(func_inputs, E_hessian1.as_vec_basic(), symbolic_cse, opt_level);

    // write another func
    parallelTransport(d1_10, t10, t1, d1_1);
    cross(t1, d1_1, d1_2);

    d2_1 = d2_10;
    cross(t2, d2_1, d2_2);

    d1_1.mul_scalar(cos(theta_1), temp1);
    d1_2.mul_scalar(sin(theta_1), temp2);
    temp1.add_matrix(temp2, m1_1);

    d1_1.mul_scalar(sin(theta_1), temp1);
    d1_2.mul_scalar(cos(theta_1), temp2);
    subtract_matrix(temp2, temp1, m1_2);

    d2_1.mul_scalar(cos(theta_2), temp1);
    d2_2.mul_scalar(sin(theta_2), temp2);
    temp1.add_matrix(temp2, m2_1);

    d2_1.mul_scalar(sin(theta_2), temp1);
    d2_2.mul_scalar(cos(theta_2), temp2);
    subtract_matrix(temp2, temp1, m2_2);

    m1_1_func_2.init(func_inputs, m1_1.as_vec_basic(), symbolic_cse, opt_level);
    m2_1_func_2.init(func_inputs, m2_1.as_vec_basic(), symbolic_cse, opt_level);
    m1_2_func_2.init(func_inputs, m1_2.as_vec_basic(), symbolic_cse, opt_level);
    m2_2_func_2.init(func_inputs, m2_2.as_vec_basic(), symbolic_cse, opt_level);
    
    // define parallel item
    m1_1.elementwise_mul_matrix(m2_1, temp1);
    parallel1 = pow(sub(add(temp1.as_vec_basic()), one), 2);    
    m1_2.elementwise_mul_matrix(m2_2, temp1);
    parallel2 = pow(sub(add(temp1.as_vec_basic()), one), 2);
    t1.elementwise_mul_matrix(t2, temp1);
    parallel3 = pow(sub(add(temp1.as_vec_basic()), one), 2);
    parallel = add(add(parallel1, parallel2), parallel3);
    
    // define the penalty energy
    x1s.add_matrix(x1e, x1);
    x1.mul_scalar(div(one, two), x1);

    x2s.add_matrix(x2e, x2);
    x2.mul_scalar(div(one, two), x2);

    subtract_matrix(x1, x2, dp);

    dp.elementwise_mul_matrix(m1_1, temp1);
    trans1 = pow(sub(add(temp1.as_vec_basic()), h1), 2);
    dp.elementwise_mul_matrix(m1_2, temp1);
    trans2 = pow(add(temp1.as_vec_basic()), 2);
    dp.elementwise_mul_matrix(t1, temp1);
    trans3 = pow(add(temp1.as_vec_basic()), 2);
    trans = add(add(trans1, trans2), trans3);

    // E = mul(div(one, two), mul(K, add(trans, parallel)));
    E = mul(div(one, two), mul(K, trans));

    DenseMatrix E_potential2{{E}};
    DenseMatrix E_gradient2(1, 14);
    DenseMatrix E_hessian2(14, 14);
    jacobian(E_potential2, nodes, E_gradient2);    
    jacobian(E_gradient2, nodes, E_hessian2);
    
    E_func_2.init(func_inputs, E_potential2.as_vec_basic(), symbolic_cse, opt_level);
    E_grad_dx_func_2.init(func_inputs, E_gradient2.as_vec_basic(), symbolic_cse, opt_level);
    E_hess_dx_func_2.init(func_inputs, E_hessian2.as_vec_basic(), symbolic_cse, opt_level);

    // write another func
    d1_1 = d1_10;
    cross(t1, d1_1, d1_2);

    d2_1 = d2_10;
    cross(t2, d2_1, d2_2);

    d1_1.mul_scalar(cos(theta_1), temp1);
    d1_2.mul_scalar(sin(theta_1), temp2);
    temp1.add_matrix(temp2, m1_1);

    d1_1.mul_scalar(sin(theta_1), temp1);
    d1_2.mul_scalar(cos(theta_1), temp2);
    subtract_matrix(temp2, temp1, m1_2);

    d2_1.mul_scalar(cos(theta_2), temp1);
    d2_2.mul_scalar(sin(theta_2), temp2);
    temp1.add_matrix(temp2, m2_1);

    d2_1.mul_scalar(sin(theta_2), temp1);
    d2_2.mul_scalar(cos(theta_2), temp2);
    subtract_matrix(temp2, temp1, m2_2);

    m1_1_func_3.init(func_inputs, m1_1.as_vec_basic(), symbolic_cse, opt_level);
    m2_1_func_3.init(func_inputs, m2_1.as_vec_basic(), symbolic_cse, opt_level);
    m1_2_func_3.init(func_inputs, m1_2.as_vec_basic(), symbolic_cse, opt_level);
    m2_2_func_3.init(func_inputs, m2_2.as_vec_basic(), symbolic_cse, opt_level);
    
    // define parallel item
    m1_1.elementwise_mul_matrix(m2_1, temp1);
    parallel1 = pow(sub(add(temp1.as_vec_basic()), one), 2);    
    m1_2.elementwise_mul_matrix(m2_2, temp1);
    parallel2 = pow(sub(add(temp1.as_vec_basic()), one), 2);
    t1.elementwise_mul_matrix(t2, temp1);
    parallel3 = pow(sub(add(temp1.as_vec_basic()), one), 2);
    parallel = add(add(parallel1, parallel2), parallel3);
    
    // define the penalty energy
    x1s.add_matrix(x1e, x1);
    x1.mul_scalar(div(one, two), x1);

    x2s.add_matrix(x2e, x2);
    x2.mul_scalar(div(one, two), x2);

    subtract_matrix(x1, x2, dp);

    dp.elementwise_mul_matrix(m1_1, temp1);
    trans1 = pow(sub(add(temp1.as_vec_basic()), h1), 2);
    dp.elementwise_mul_matrix(m1_2, temp1);
    trans2 = pow(add(temp1.as_vec_basic()), 2);
    dp.elementwise_mul_matrix(t1, temp1);
    trans3 = pow(add(temp1.as_vec_basic()), 2);
    trans = add(add(trans1, trans2), trans3);
    
    E = mul(div(one, two), mul(K, trans));
    // E = mul(div(one, two), mul(K, add(trans, parallel)));

    DenseMatrix E_potential3{{E}};
    DenseMatrix E_gradient3(1, 14);
    DenseMatrix E_hessian3(14, 14);
    jacobian(E_potential3, nodes, E_gradient3);    
    jacobian(E_gradient3, nodes, E_hessian3);

    E_func_3.init(func_inputs, E_potential3.as_vec_basic(), symbolic_cse, opt_level);
    E_grad_dx_func_3.init(func_inputs, E_gradient3.as_vec_basic(), symbolic_cse, opt_level);
    E_hess_dx_func_3.init(func_inputs, E_hessian3.as_vec_basic(), symbolic_cse, opt_level);
    
    


}


   



// void symbolicEquations::generatePenaltyFunctions(){
//     DenseMatrix x1s({x1s_x, x1s_y, x1s_z});
//     DenseMatrix x1e({x1e_x, x1e_y, x1e_z});
//     DenseMatrix x2s({x2s_x, x2s_y, x2s_z});
//     DenseMatrix x2e({x2e_x, x2e_y, x2e_z});

//     DenseMatrix d1_10({d1_1x0, d1_1y0, d1_1z0});
//     DenseMatrix d1_20({d1_2x0, d1_2y0, d1_2z0});
//     DenseMatrix d2_10({d2_1x0, d2_1y0, d2_1z0});
//     DenseMatrix d2_20({d2_2x0, d2_2y0, d2_2z0});



//     int num_rows = x1s.nrows();
//     int num_cols = x1s.ncols();
    
//     DenseMatrix temp1(3, 1);
//     DenseMatrix temp2(3, 1);

//     DenseMatrix m1_1(3, 1);
//     DenseMatrix m2_1(3, 1);

//     d1_10.mul_scalar(cos(theta_1), temp1);
//     d1_20.mul_scalar(sin(theta_1), temp2);
//     temp1.add_matrix(temp2, m1_1);


//     DenseMatrix x1(3, 1);
//     x1s.add_matrix(x1e, x1);
//     x1.mul_scalar(div(one, two), x1);

//     m1_1.mul_scalar(div(h1, two), m1_1);
//     subtract_matrix(x1, m1_1, x1);
    

//     DenseMatrix x2(3, 1);
//     d2_10.mul_scalar(cos(theta_2), temp1);
//     d2_20.mul_scalar(sin(theta_2), temp2);
//     temp1.add_matrix(temp2, m2_1);
//     x2s.add_matrix(x2e, x2);
//     x2.mul_scalar(div(one, two), x2);

//     m2_1.mul_scalar(div(h1, two), m2_1);
//     x2.add_matrix(m2_1, x2);
    

//     DenseMatrix x_rel(3, 1);
//     subtract_matrix(x2, x1, x_rel);

//     RCP<const Basic> x_rel_norm;
//     get_norm(x_rel, x_rel_norm);
    
//     RCP<const Basic> E = mul(div(one, two), mul(K, pow(x_rel_norm, 2)));

//     DenseMatrix E_potential{{E}};

//     DenseMatrix E_gradient(1, 14);
//     DenseMatrix E_hessian(14, 14);

//     vec_basic nodes_vec {x1s_x, x1s_y, x1s_z, theta_1,
//                          x1e_x, x1e_y, x1e_z,
//                          x2s_x, x2s_y, x2s_z, theta_2,
//                          x2e_x, x2e_y, x2e_z};

//     DenseMatrix nodes {nodes_vec};
//     jacobian(E_potential, nodes, E_gradient);    
//     jacobian(E_gradient, nodes, E_hessian);

//     vec_basic func_inputs{x1s_x, x1s_y, x1s_z, theta_1,
//                           x1e_x, x1e_y, x1e_z,
//                           x2s_x, x2s_y, x2s_z, theta_2,
//                           x2e_x, x2e_y, x2e_z,
//                           d1_1x0, d1_1y0, d1_1z0,
//                           d1_2x0, d1_2y0, d1_2z0,
//                           d2_1x0, d2_1y0, d2_1z0,
//                           d2_2x0, d2_2y0, d2_2z0};
    
//     func_inputs.push_back(K);
//     func_inputs.push_back(h1);

//     // E_func.init(func_inputs, x1.as_vec_basic(), symbolic_cse, opt_level);
//     E_func.init(func_inputs, E_potential.as_vec_basic(), symbolic_cse, opt_level);
//     E_grad_dx_func.init(func_inputs, E_gradient.as_vec_basic(), symbolic_cse, opt_level);
//     E_hess_dx_func.init(func_inputs, E_hessian.as_vec_basic(), symbolic_cse, opt_level);


//     // std::cout << E_potential.nrows()<<" "<<E_potential.ncols() << std::endl;


//     vec_basic frames_vec {d1_1x0, d1_1y0, d1_1z0,
//                           d1_2x0, d1_2y0, d1_2z0,
//                           d2_1x0, d2_1y0, d2_1z0,
//                           d2_2x0, d2_2y0, d2_2z0};
    
//     DenseMatrix frames {frames_vec};
    
//     DenseMatrix E_frame_grad(1, 12);
//     DenseMatrix E_frame_hess(12, 12);

//     jacobian(E_potential, frames, E_frame_grad);    
//     jacobian(E_frame_grad, frames, E_frame_hess);

//     E_grad_dd_func.init(func_inputs, E_frame_grad.as_vec_basic(), symbolic_cse, opt_level);
//     E_hess_dd_func.init(func_inputs, E_frame_hess.as_vec_basic(), symbolic_cse, opt_level);
    
//     m1_1_func.init(func_inputs, m1_1.as_vec_basic(), symbolic_cse, opt_level);
//     m2_1_func.init(func_inputs, m1_1.as_vec_basic(), symbolic_cse, opt_level);

// }
