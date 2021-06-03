//
// Created by kazem on 5/31/21.
//


#include <vector>
#include <thread>

#include "SympilerSolver.hpp"

namespace IPC{


 template <typename vectorTypeI, typename vectorTypeS>
 SympilerSolver<vectorTypeI, vectorTypeS>::SympilerSolver(void)
 {
  A = new sym_lib::parsy::CSC;
  sym_chol = NULL;
 }

 template <typename vectorTypeI, typename vectorTypeS>
 SympilerSolver<vectorTypeI, vectorTypeS>::~SympilerSolver(void)
 {
  delete A;
  delete sym_chol;
 }

 template <typename vectorTypeI, typename vectorTypeS>
 void SympilerSolver<vectorTypeI, vectorTypeS>::set_pattern(const std::vector<std::set<int>>& vNeighbor,
                                                           const std::set<int>& fixedVert)
 {
  Base::set_pattern(vNeighbor, fixedVert);
  Base::ia.array() -= 1;
  Base::ja.array() -= 1; // CHOLMOD's index starts from 0
  A->i = Base::ja.data();
  A->p = Base::ia.data();
  A->x = Base::a.data();
  A->nzmax = Base::ja.size();
  A->ncol = A->nrow = Base::numRows;
  A->packed = 1;
  A->stype = -1; // lower triangular
 }

 template <typename vectorTypeI, typename vectorTypeS>
 void SympilerSolver<vectorTypeI, vectorTypeS>::set_pattern(const Eigen::SparseMatrix<double>& mtr)
 {
  Base::set_pattern(mtr);

   Ax = A->x;
   Ap = A->p;
   Ai = A->i;
   // -1: upper right part will be ignored during computation

   A->i = Base::ja.data();
   A->p = Base::ia.data();
   A->x = Base::a.data();
  A->nzmax = Base::ja.size();
  A->ncol = A->nrow = Base::numRows;
  A->packed = 1;
  A->stype = -1; // lower triangular
 }

 template <typename vectorTypeI, typename vectorTypeS>
 void SympilerSolver<vectorTypeI, vectorTypeS>::load(const char* filePath, Eigen::VectorXd& rhs)
 {
  Base::load(filePath, rhs);

  Base::ia.array() -= 1;
  Base::ja.array() -= 1; // CHOLMOD's index starts from 0
  A->i = Base::ja.data();
  A->p = Base::ia.data();
  A->x = Base::a.data();
  A->nzmax = Base::ja.size();
  A->ncol = A->nrow = Base::numRows;
  A->packed = 1;
  A->stype = -1; // lower triangular
 }

 template <typename vectorTypeI, typename vectorTypeS>
 void SympilerSolver<vectorTypeI, vectorTypeS>::analyze_pattern(void)
 {
  // std::cout << getCurrentRSS() << std::endl;
  //cholmod_free_factor(&L, &cm);
  //L = cholmod_analyze(A, &cm);
  delete sym_chol;
  const auto processor_count = std::thread::hardware_concurrency();
  sym_chol = sympiler::sympiler_chol_symbolic(A, processor_count);
  //std::cout<<"\n*********************\n";
 }

 template <typename vectorTypeI, typename vectorTypeS>
 bool SympilerSolver<vectorTypeI, vectorTypeS>::factorize(void)
 {
  sym_chol->numerical_factorization(A);
  //cholmod_factorize(A, L, &cm);
  // std::cout << getCurrentRSS() << std::endl;
  // exit(0);
 }

 template <typename vectorTypeI, typename vectorTypeS>
 void SympilerSolver<vectorTypeI, vectorTypeS>::solve(Eigen::VectorXd& rhs,
                                                     Eigen::VectorXd& result)
 {
  //TODO: directly point to rhs?
  //b->x = rhs.data();
  //cholmod_dense* x;
  //x = cholmod_solve(CHOLMOD_A, L, b, &cm);
  result.conservativeResize(rhs.size());
  //memcpy(result.data(), x->x, result.size() * sizeof(result[0]));
  //cholmod_free_dense(&x, &cm);
  auto *sol_dbl = sym_chol->solve_only(rhs.data(),1);
  //for (int i = 0; i < rhs.size(); ++i) {
 //  std::cout<<sol_dbl[i]<<",";
 // }
  //std::cout<<"\n";
  //sym_chol->compute_norms(rhs.data());
  //std::cout<<" ====> "<<sym_chol->res_l1<<"\n";
  result.conservativeResize(rhs.size());
  memcpy(result.data(), sol_dbl, rhs.size() * sizeof(double));
  //result = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> >(
  //  sol_dbl,A->ncol,1);
 }

 template <typename vectorTypeI, typename vectorTypeS>
 void SympilerSolver<vectorTypeI, vectorTypeS>::multiply(const Eigen::VectorXd& x,
                                                        Eigen::VectorXd& Ax)
 {
  assert(x.size() == Base::numRows);

 }

 template <typename vectorTypeI, typename vectorTypeS>
 void SympilerSolver<vectorTypeI, vectorTypeS>::outputFactorization(const std::string& filePath)
 {
  //cholmod_sparse* spm = cholmod_factor_to_sparse(L, &cm);

  FILE* out = fopen(filePath.c_str(), "w");
  assert(out);

  //cholmod_write_sparse(out, spm, NULL, "", &cm);

  fclose(out);
 }

 template class SympilerSolver<Eigen::VectorXi, Eigen::VectorXd>;

}