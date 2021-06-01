//
// Created by kazem on 5/31/21.
//

#ifndef IPC_SYMPILERSOLVER_HPP
#define IPC_SYMPILERSOLVER_HPP


#include "LinSysSolver.hpp"

#include <sympiler_cholesky.h>

namespace IPC{

 template <typename vectorTypeI, typename vectorTypeS>
 class SympilerSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
  typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:

  sym_lib::parsy::CSC *A;
  sym_lib::parsy::SolverSettings *sym_chol;
  int *Ai, *Ap;
  double *Ax, *bx, *solutionx, *x_cdx, *y_cdx;

 public:
  SympilerSolver(void);
  ~SympilerSolver(void);

  void set_pattern(const std::vector<std::set<int>>& vNeighbor,
                   const std::set<int>& fixedVert);
  void set_pattern(const Eigen::SparseMatrix<double>& mtr); //NOTE: mtr must be SPD
  void load(const char* filePath, Eigen::VectorXd& rhs);

  void analyze_pattern(void);

  bool factorize(void);

  void solve(Eigen::VectorXd& rhs,
             Eigen::VectorXd& result);

  virtual void multiply(const Eigen::VectorXd& x,
                        Eigen::VectorXd& Ax);

  virtual void outputFactorization(const std::string& filePath);
 };

}

#endif //IPC_SYMPILERSOLVER_HPP
