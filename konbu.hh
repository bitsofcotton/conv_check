/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/* BSD 3-Clause License:
 * Copyright (c) 2013 - 2021, kazunobu watatsu.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(_LINEAR_OPT_)

using std::cerr;
using std::flush;

template <typename T> SimpleVector<T> inner(const SimpleMatrix<T>& A, const SimpleVector<T>& bl, const SimpleVector<T>& bu) {
#if defined(_FLOAT_BITS_)
  static const auto epsilon(T(1) >> int64_t(mybits / 2));
#else
  static const auto epsilon(sqrt(std::numeric_limits<T>::epsilon()));
#endif
  assert(A.rows() == bl.size() && A.rows() == bu.size() &&
         0 < A.cols() && 0 < A.rows());
  // cout << A << endl << b << endl << c.transpose() << endl;
  cerr << " (" << A.rows() << ", " << A.cols() << ")";
  // bu - bb == A, bl - bb == - A <=> bu - bl == 2 A
  const auto bb(bu - (bu - bl) / T(2));
  const auto upper(bu - bb);
  SimpleMatrix<T> AA(A.rows() * 2 - 1 + (A.cols() + 1) * 2 - 1, A.cols() + 1);
  SimpleVector<T> one(AA.rows());
  for(int i = 0; i < A.rows(); i ++) {
    for(int j = 0; j < A.cols(); j ++)
      AA(i, j) = A(i, j);
    AA(i, A.cols()) = - bb[i];
    if(upper.dot(upper) != T(0)) {
      assert(upper[i] != T(0));
      AA.row(i) /= upper[i];
    }
    one[i] = T(1);
    if(i < A.rows() - 1) {
      AA.row(i + A.rows()) = - AA.row(i);
      one[i + A.rows()] = T(1);
    }
    assert(isfinite(AA.row(i).dot(AA.row(i))));
  }
  // [A, - bb] [x t] <= [1 ... 1].
  for(int i = A.rows() * 2 - 1; i < AA.rows(); i ++) {
    for(int j = 0; j < AA.cols(); j ++)
      AA(i, j) = (i - (A.rows() * 2 - 1)) / 2 == j
        ? T((i - (A.rows() * 2 - 1)) & 1 ? 1 : - 1)
        : T(0);
    one[i] = T(1);
  }
  SimpleMatrix<T> Pt(AA.cols(), AA.rows());
  for(int i = 0; i < Pt.rows(); i ++)
    for(int j = 0; j < Pt.cols(); j ++)
      Pt(i, j) = T(0);
  for(int i = 0; i < AA.cols(); i ++) {
    const auto Atrowi(AA.col(i));
    const auto work(Atrowi - Pt.projectionPt(Atrowi));
    // generally, assert norm > error is needed.
    // in this case, not.
    Pt.row(i) = work / sqrt(work.dot(work));
  }
  cerr << "Q" << flush;
  const auto R(Pt * AA);
  cerr << "R" << flush;
  for(int n_fixed = 0; n_fixed < Pt.rows() - 1; n_fixed ++) {
    const auto on(Pt.projectionPt(- one));
    auto fidx(0);
    for( ; fidx < on.size() && on[fidx] <= T(0); fidx ++) ;
    for(int i = fidx + 1; i < on.size(); i ++)
      if(T(0) < on[i] && on[i] < on[fidx]) fidx = i;
    if(! n_fixed) fidx = on.size() - 1;
    if(on.size() <= fidx || on[fidx] <= T(0)) break;
    
    // O(mn^2) over all in this function.
    const auto orth(Pt.col(fidx));
    const auto norm2orth(orth.dot(orth) + (n_fixed ? epsilon : T(1)));
    if(norm2orth <= T(0)) break;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * (Pt.col(j).dot(orth) + (n_fixed ? epsilon : T(1))) / norm2orth);
  }
  cerr << "G" << flush;
#if defined(_WITHOUT_EIGEN_)
  auto rvec(- R.solve(Pt * one));
#else
  auto rvec(- R.inverse() * (Pt * one));
#endif
  cerr << "I" << flush;
  SimpleVector<T> rrvec(rvec.size() - 1);
  // | [A, - bb == - upper] [x t] | <= epsilon 1.
  for(int i = 0; i < rrvec.size(); i ++)
    rrvec[i] = rvec[i] / rvec[rvec.size() - 1];
  return rrvec;
}

#define _LINEAR_OPT_
#endif

