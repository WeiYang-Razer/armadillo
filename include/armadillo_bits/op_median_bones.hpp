// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup op_median
//! @{


template<typename T>
struct arma_cx_median_packet
  {
  T     val;
  uword index;
  
  inline arma_cx_median_packet() : val(0), index(0) {}
  
  inline explicit arma_cx_median_packet(const T in_val, const uword in_index) : val(in_val), index(in_index) {}
  inline          arma_cx_median_packet(const arma_cx_median_packet& in     ) : val(in.val), index(in.index) {}
  
  inline void operator=(const arma_cx_median_packet& in) { if(this == &in) { return; } val = in.val; index = in.index; }
  };



template<typename T>
arma_inline
bool
operator< (const arma_cx_median_packet<T>& A, const arma_cx_median_packet<T>& B)
  {
  return (A.val < B.val);
  }



class op_median
  : public traits_op_xvec
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_median>& expr);
  
  template<typename eT>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim, const typename arma_not_cx<eT>::result* junk = nullptr);
  
  template<typename eT>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim, const typename arma_cx_only<eT>::result* junk = nullptr);
  
  //
  //
  
  template<typename T1>
  inline static typename T1::elem_type median_vec(const T1& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = nullptr);
  
  template<typename T1>
  inline static typename T1::elem_type median_vec(const T1& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr);
  
  //
  //
  
  template<typename eT>
  inline static eT direct_median(std::vector<eT>& X);
  
  template<typename T>
  inline static void direct_cx_median_index(uword& out_index1, uword& out_index2, std::vector< arma_cx_median_packet<T> >& X);
  };



class op_median_omit
  : public traits_op_xvec
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_median_omit>& expr);
  
  template<typename eT, typename functor>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim, functor is_omitted, const typename arma_not_cx<eT>::result* junk = nullptr);
  
  template<typename eT, typename functor>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim, functor is_omitted, const typename arma_cx_only<eT>::result* junk = nullptr);
  
  template<typename T1, int omit_mode>
  inline static typename T1::elem_type median_vec(const T1& X, const elem_opts::omit_indicator<omit_mode>&, const typename arma_not_cx<typename T1::elem_type>::result* junk = nullptr);
  
  template<typename T1, int omit_mode>
  inline static typename T1::elem_type median_vec(const T1& X, const elem_opts::omit_indicator<omit_mode>&, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr);
  };



//! @}
