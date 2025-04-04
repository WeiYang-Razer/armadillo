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


//! \addtogroup fn_balance
//! @{



template<typename T1>
inline
bool
balance(Col<typename T1::pod_type>& D, Col<uword>& P, Mat<typename T1::elem_type>& B, const Base<typename T1::elem_type,T1>& A, const char* operation = "both")
  {
  arma_debug_sigprint();
  
  // TODO: check for aliasing between D and B
  
  const char sig = (operation != nullptr) ? operation[0] : char(0);
  
  if( (sig != 'b') && (sig != 'p') && (sig != 's') )  { arma_stop_logic_error("balance(): unsupported operation"); }
  
  const bool do_perm  = (sig == 'b') || (sig == 'p');
  const bool do_scale = (sig == 'b') || (sig == 's');
  
  const bool calc_DP = true;
  
  B = A.get_ref();
  
  if(B.is_square() == false)
    {
    arma_stop_logic_error("balance(): given matrix must be square sized");
    
    return false;
    }
  
  const bool status = auxlib::balance(D, P, B, calc_DP, do_perm, do_scale);
  
  if(status == false)
    {
    // TODO: emit warning
    D.soft_reset();
    P.soft_reset();
    B.soft_reset();
    }
  
  return status;
  }



template<typename T1>
inline
bool
balance(Mat<typename T1::elem_type>& B, const Base<typename T1::elem_type,T1>& A, const char* operation = "both")
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const char sig = (operation != nullptr) ? operation[0] : char(0);
  
  if( (sig != 'b') && (sig != 'p') && (sig != 's') )  { arma_stop_logic_error("balance(): unsupported operation"); }
  
  const bool do_perm  = (sig == 'b') || (sig == 'p');
  const bool do_scale = (sig == 'b') || (sig == 's');
  
  const bool calc_DP = false;
  
  Col<T>     D;
  Col<uword> P;
  
  B = A.get_ref();
  
  if(B.is_square() == false)
    {
    arma_stop_logic_error("balance(): given matrix must be square sized");
    
    return false;
    }
  
  const bool status = auxlib::balance(D, P, B, calc_DP, do_perm, do_scale);
  
  if(status == false)
     {
     // TODO: emit warning
     B.soft_reset();
     }
  
  return status;
  }



template<typename T1>
inline
Mat<typename T1::elem_type>
balance(const Base<typename T1::elem_type,T1>& A, const char* operation = "both")
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> B;
  
  const bool status = balance(B, A, operation);
  
  if(status == false)  { arma_stop_runtime_error("balance(): transformation failed"); }
  
  return B;
  }



//! @}
