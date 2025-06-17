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


//! \addtogroup op_omit
//! @{



template<typename T1>
inline
void
op_omit::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_omit>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword omit_mode = in.aux_uword_a;
  
  auto is_omitted_1 = [](const eT& x) -> bool { return arma_isnan(x);       };
  auto is_omitted_2 = [](const eT& x) -> bool { return arma_isnonfinite(x); };
  
  if(omit_mode == 1)  { op_omit::apply(out, in.m, is_omitted_1); }
  if(omit_mode == 2)  { op_omit::apply(out, in.m, is_omitted_2); }
  }



template<typename T1, typename functor>
inline
void
op_omit::apply(Mat<typename T1::elem_type>& out, const T1& X, functor is_omitted)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Mat<T1>::value || is_subview_col<T1>::value || is_Mat<typename Proxy<T1>::stored_type>::value || Proxy<T1>::use_mp)
    {
    const quasi_unwrap<T1> U(X);
    
    const eT*   X_mem = U.M.memptr();
    const uword N     = U.M.n_elem;
    
    Mat<eT> Y(N, 1, arma_nozeros_indicator());
    
    eT* Y_mem = Y.memptr();
    
    uword count = 0;
    
    for(uword i=0; i < N; ++i)
      {
      const eT val = X_mem[i];
      
      if(is_omitted(val) == false)  { Y_mem[count] = val; ++count; }
      }
    
    out.steal_mem_col(Y, count);
    }
  else
    {
    const Proxy<T1> P(X);
    
    const uword N = P.get_n_elem();
    
    Mat<eT> Y(N, 1, arma_nozeros_indicator());
    
    eT* Y_mem = Y.memptr();
    
    uword count = 0;
    
    if(Proxy<T1>::use_at == false)
      {
      const typename Proxy<T1>::ea_type Pea = P.get_ea();
      
      for(uword i=0; i < N; ++i)
        {
        const eT val = Pea[i];
        
        if(is_omitted(val) == false)  { Y_mem[count] = val; ++count; }
        }
      }
    else
      {
      const uword n_rows = P.get_n_rows();
      const uword n_cols = P.get_n_cols();
      
      for(uword c=0; c < n_cols; ++c)
      for(uword r=0; r < n_rows; ++r)
        {
        const eT val = P.at(r,c);
        
        if(is_omitted(val) == false)  { Y_mem[count] = val; ++count; }
        }
      }
      
    out.steal_mem_col(Y, count);
    }
  }



//! @}
