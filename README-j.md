# Konbu Check
複数の線形制約条件から内点を一つとってくるプログラムです。O(mn^2) 回の四則演算及び初等関数演算が必要になります。
また、もっと早く計算するためには、Eigen ライブラリが必要です。さらに精度よく計算するためには、mpfr++ ライブラリと付随するライブラリ、
または、QD ライブラリが必要です。  
現在開発休止中です。

# Tips
'err_error' または 'intercept' に続く stderr に出力される数値は問題に依存していて、このプログラムがその精度内で
解くことが難しいかどうかの指標になります。値が 0 よりも極端に大きい場合 (特に >= 1 の時)には、得られる値は非常に訝しいものとなります。  
また、元の問題の実行可能領域が極端に薄い場合には、とってきた値が意味をなさないことがあります。
その場合には、パラメータ threshold_loop を少し広げた範囲でもう一度実行してみてください。

# バグ
パラメタか精度が問題に対して悪い値の時に、とってくる値がおかしくなることがあります。
よくスケールされた問題に対してはこれは起こりにくいです。  
また、最適化に使用する際には、パラメータ調整内部で、最適解の変数に対する大きさの比率の上限をある程度規定しています。
これにより、あまり最適でない解が得られることがあります。

# 調整可能なパラメタ
konbu.hh 内の LP<T>::LP() は調整可能です。
* threshold_feas   : QR 分解の誤差です。
* threshold_p0     : ループ時の P 行列の誤差です。
* threshold_loop   : [-b+1&epsilon;,P][t,x]&leq;b, &epsilon;
* threshold_inner  : ループ時に内点かどうか判断する誤差です。
* largest_intercept: Ax&leq;b を Ax&geq;b を除いてとってくる際に使う超立体の大きさです。
* largest_opt      : 最適解の最大の大きさの比率を仮定します。
* n_opt_steps      : 最適解をとってくる際のループ回数です。

# 文脈
特許出願の JP2014089683 . 

# 開発状態
g++ または clang++ が OpenACC 経由で iris チップをサポートすることを待っています。  

# インストール方法
Makefile を使用するライブラリが適切に通るように変更してコンパイルしてください。
オプションには、-DACC_GMP=$bits オプション及び、-DACC_QD_QDOUBLE オプション、または -DACC_DOUBLE オプションなどのうち
一つを指定してください。また、-DWITHOUT_EIGEN オプションは非常に遅いです。

# デモ
https://services.limpid-intensity.info/konbu.php にあります。

# 証明
Ax&lt;=b

[-b,P][t,x']&lt;=0,
P is part of orthogonal matrix.

[-b'+1&epsilon;,P][t,x'+b'']&lt;=1&epsilon;,
b' is orthogonal to P.

After loop we get :
P'[Q[t,x'',0]]&lt;=0

# 使い方
    #include "konbu_init.h"
    ...
    // if you use with mpfr, num_t::set_default_proc(BITS); is needed.
    // if you use with QD,   unsigned int old_cw; fpu_fix_start(&old_cw); is needed.
    ...
    int m; // number of rows;
    int n; // number of columns;
    Mat A(m, n);
    ...
    Vec result;
    bool fix_partial[A.rows()];
    LP<num_t> lp;
    bool feas = lp.inner(fix_partial, result, A, b);
    // if you use with QD,    fpu_fix_end(&old_cw); is needed.

# 少しの情報
パターンマッチ(補空間の&infin;-ノルム最小化)にも使えます。
32 ビット機で使う場合には、メモリの関係でおおよそ 8k x 4k のサイズまでしか扱えません。  
もしも、解くべき問題の変数の最大最小値が [-&pi;,&pi;] 範囲内にあるなど、安定している場合には、
largest_intercept パラメータをそのように変えることで、(低い精度でも)かなり安定した結果を得ることができます。

# その他のダウンロードサイト
* https://ja.osdn.net/projects/conv-check/
* https://www.sourceforge.net/projects/convcheck/
* https://konbu.sakura.ne.jp/files/konbu_check-1.01-release.tar.gz
* https://files.limpid-intensity.info/konbu_check-1.01-release.tar.gz
