#include <cstdio>
#include <cstdlib>
#include "konbu_init.h"

#define SZ_NUM_BUF 500

int main(int argc, char* argv[])
{
  FILE* in;
  in = fopen(argv[1], "ro");
  if(!in) return -1;

#if defined(ACC_GMP)
  num_t::set_default_prec(mpfr::digits2bits(ACC_GMP));
#elif defined(ACC_QD_DDOUBLE) || defined(ACC_QD_QDOUBLE)
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
#endif
  
  cout << "File: " << argv[1] << endl;
  cout.precision(MANT_DIG);
  cerr.precision(MANT_DIG);
  
  int m, n;
  
  fscanf(in, "%d", &m);
  fscanf(in, "%d", &n);
  
  Mat A(n, m);
  Vec b(n);
  char numbuf[SZ_NUM_BUF];
  
  for(int i = 0; i < n; i ++)
    for(int j = 0; j < m; j ++){
      bool ws_loop = true;
      while(ws_loop) {
        char c = fgetc(in);
        switch(c) {
          case ' ': case '\t': case '\n':
            break;
          default:
            ws_loop = false;
            fseek(in, -1, SEEK_CUR);
        }
      }
      bool num_loop = true;
      for(int k = 0; k < SZ_NUM_BUF && num_loop;) {
        char cc = fgetc(in);
        switch(cc) {
          case '0': case '1': case '2': case '3': case '4':
          case '5': case '6': case '7': case '8': case '9':
          case '.': case '-': case 'e': case '+':
            numbuf[k ++] = cc;
            break;
          default:
            numbuf[k ++] = 0;
            fseek(in, -1, SEEK_CUR);
            num_loop  = false;
        }
      }
#if defined(ACC_GMP)
      A(i, j) = std::string(numbuf);
#else
      std::stringstream(std::string(numbuf)) >> A(i, j);
#endif
    }
  for(int i = 0; i < n; i ++) {
    bool ws_loop = true;
    while(ws_loop) {
      char c = fgetc(in);
      switch(c) {
        case ' ': case '\t': case '\n':
          break;
        default:
          ws_loop = false;
          fseek(in, -1, SEEK_CUR);
      }
    }
    bool num_loop = true;
    for(int j = 0; j < SZ_NUM_BUF && num_loop;) {
      char cc = fgetc(in);
      switch(cc) {
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
        case '.': case '-': case 'e': case '+':
          numbuf[j ++] = cc;
          break;
        default:
          numbuf[j ++] = 0;
          fseek(in, -1, SEEK_CUR);
          num_loop  = false;
      }
    }
#if defined(ACC_GMP)
    b[i] = std::string(numbuf);
#else
    std::stringstream(std::string(numbuf)) >> b[i];
#endif
  }
  fclose(in);
  
  cerr << " PARSE";
  fflush(stdout);
  
  Vec  result;
  bool* fix_partial = new bool[A.rows()];
  LP<num_t> lp;
  bool feas = lp.inner(fix_partial, result, A, b);
  
  for(int i = 0; i < result.size(); i ++)
    cout << result[i] << ", ";
  cout << endl;
  
  int n_fixed = 0;
  for(int i = 0; i < A.rows(); i ++)
    if(fix_partial[i]) n_fixed ++;
  cout << "fix_partial ( " << n_fixed << " / " << A.cols() << " ) : ";
  for(int i = 0; i < A.rows(); i ++)
    cout << (const char*)(fix_partial[i] ? "1" : "0");
  cout << endl;
  
  delete[] fix_partial;
  
#if defined(ACC_QD_DDOUBLE) || defined(ACC_QD_QDOUBLE)
  fpu_fix_end(&old_cw);
#endif
  
  return 0;
}

