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
  static const auto epsilon(T(1) >> int64_t(mybits - 1));
#else
  static const auto epsilon(std::numeric_limits<T>::epsilon());
#endif
  static const auto ee(pow(epsilon, T(1) / T(8)));
  assert(A.rows() == bl.size() && A.rows() == bu.size() &&
         0 < A.cols() && 0 < A.rows());
  // cout << A << endl << b << endl << c.transpose() << endl;
  cerr << " (" << A.rows() << ", " << A.cols() << ")";
  // bu - bb == A, bl - bb == - A <=> bu - bl == 2 A
  const auto bb((bu + bl) / T(2));
  const auto upper(bu - bb);
  SimpleMatrix<T> AA(A.rows() * 2 - 1 + A.cols() * 2, A.cols() + 1);
  SimpleVector<T> one(AA.rows());
  std::vector<std::pair<T, int> > fidx0;
  fidx0.reserve(A.rows());
  for(int i = 0; i < A.rows(); i ++) {
    for(int j = 0; j < A.cols(); j ++)
      AA(i, j) = A(i, j);
    AA(i, A.cols()) = - bb[i];
    one[i]          = T(1);
    assert(T(0) <= upper[i]);
    if(upper[i] == T(0)) {
      const auto n2(AA.row(i).dot(AA.row(i)));
      if(n2 != T(0)) {
        fidx0.emplace_back(std::make_pair(- T(1), i));
        AA.row(i) /= sqrt(n2);
      } else
        AA.row(i) *= n2;
    } else
      AA.row(i) /= upper[i];
    const auto A2(AA.row(i).dot(AA.row(i)));
    assert(isfinite(A2) && ! isnan(A2));
    if(A.rows() - 1 <= i) break;
    AA.row(i + A.rows()) = - AA.row(i);
    one[i + A.rows()] = T(1);
  }
  // this is tricky but one of them is fixed if it is needed.
  for(int i = A.rows() * 2 - 1; i < AA.rows(); i ++) {
    const auto ii(i - (A.rows() * 2 - 1));
    for(int j = 0; j < A.cols(); j ++)
      AA(i, j) = T(j == ii / 2 ? (ii & 1 ? - 1 : 1) : 0);
    AA(i, A.cols()) = T(ii & 1 ? 1 : - 1);
    AA.row(i) /= ee;
    one[i] = T(1);
  }
  // N.B. we now have |[A -bb] [x t]| <= 1 condition.
  // N.B. there's no difference |[A - bb] [x t]|^2 <= 1 condition in this.
  //      but not with mixed condition.
  SimpleMatrix<T> Pt(AA.cols(), AA.rows());
  for(int i = 0; i < Pt.rows(); i ++)
    for(int j = 0; j < Pt.cols(); j ++)
      Pt(i, j) = T(0);
  std::vector<int> residue;
  for(int i = 0; i < Pt.rows(); i ++) {
    const auto Atrowi(AA.col(i));
    const auto work(Atrowi - Pt.projectionPt(Atrowi));
    const auto n2(work.dot(work));
    if(n2 <= epsilon) {
      residue.emplace_back(i);
      continue;
    }
    Pt.row(i) = work / sqrt(n2);
  }
  int ii(0);
  for(int j = 0; j < Pt.cols() && ii < residue.size(); j ++) {
    SimpleVector<T> ek(Pt.cols());
    for(int k = 0; k < Pt.cols(); k ++)
      ek[k] = j == k ? T(1) : T(0);
    ek -= Pt.projectionPt(ek);
    const auto n2(ek.dot(ek));
    if(n2 <= epsilon) continue;
    Pt.row(residue[ii ++]) = ek / sqrt(n2);
  }
  assert(residue.size() <= ii);
  cerr << "Q" << flush;
  const auto R(Pt * AA);
  cerr << "R" << flush;
  // we now have: Q [R [x t] ] <= {0, 1}^m cond.
  const auto on(Pt.projectionPt(one));
        auto fidx(fidx0);
  fidx.reserve(fidx.size() + on.size());
  for(int i = 0; i < on.size(); i ++)
    if(isfinite(on[i]) && ! isnan(on[i]))
      fidx.emplace_back(std::make_pair(abs(on[i]), i));
  if(fidx.size())
    std::sort(fidx.begin(), fidx.end());
  // sort by: |<Q^t(1), q_k>|, we subject to minimize each, to do this,
  //   maximize minimum q_k orthogonality.
  for(int n_fixed = 0, idx = 0; n_fixed < Pt.rows() - 1 && idx < fidx.size(); n_fixed ++, idx ++) {
    const auto& iidx(fidx[idx].second);
    const auto  orth(Pt.col(iidx));
    const auto  n2(orth.dot(orth));
    if(n2 <= epsilon) {
      n_fixed --;
      continue;
    }
    // N.B. O(mn) can be writed into O(lg m + lg n) in many core cond.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / n2);
  }
  // N.B. we now have: |Q x'| <= 1 * epsilon condition.
        auto rrvec(Pt * one);
  const auto r2(rrvec.dot(rrvec));
  cerr << "I(" << r2 << ", " << rrvec[rrvec.size() - 1] << ")" << flush;
  SimpleVector<T> rvec(rrvec.size() - 1);
  if(! (isfinite(r2) && ! isnan(r2))) {
    for(int i = 0; i < rvec.size(); i ++)
      rvec[i] = T(0);
    return rvec;
  }
  rrvec = R.solve(rrvec);
  //      <=> |AA rvec| = |[A -bb] [x t]| <= 1 * epsilon.
  //      <=> - 1 * |epsilon / t| <= Ax - bb <= 1 * |epsilon / t|.
  //      if we have eps := |epsilon / t| <= 1 condition, it's ok.
  for(int i = 0; i < rvec.size(); i ++)
    rvec[i] = rrvec[i];
  const auto errb(A * rvec - bb);
        T    M(0);
  for(int i = 0; i < errb.size(); i ++)
    if(upper[i] != T(0))
      M = max(M, abs(errb[i] / upper[i]));
  if(isfinite(M) && ! isnan(M) && M != T(0))
    return rvec /= M;
  return rvec;
}

#define _LINEAR_OPT_
#endif

