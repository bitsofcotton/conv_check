# Konbu Check
複数の線形制約条件から内点を一つとってくるプログラムです。m 制約式 n 変数の問題に対して O(mn^2) 回の四則演算及び O(mn) 回の平方根演算が必要になります。
また、もっと速く計算するためには、Eigen ライブラリが必要です。さらに精度よく計算するためには(ほとんどの場合には必要ですが)、
mpfr++ ライブラリと付随するライブラリ、または、QD ライブラリが必要です。  

# 文脈
特許出願の JP2014089683 . 

# インストール方法
Makefile を使用するライブラリが適切に通るように変更してコンパイルしてください。
オプションには、-DACC_GMP=$bits オプション及び、-DACC_QD_QDOUBLE オプション、または -DACC_DOUBLE オプションなどのうち
一つを指定してください。また、-D_WITHOUT_EIGEN_ または -DACC_NO_FLOAT オプションは非常に遅いです。

# 証明
Ax&lt;=b

[-b,P][t,x']&lt;=0,
P is part of orthogonal matrix.

[-b'+1&epsilon;,P][t,x'+b'']&lt;=1&epsilon;,
b' is orthogonal to P.

P'[Q[t,x'',0]]&lt;=0
P'の直交性はループ中に保存されます。

0&lt;a となるベクトルを両辺に加え、左辺の変動を無視することによって
||a||/||x||-&lt; 0 ととれる場合には良い結果を返します。
この変換は実行可能領域を変更します。
また、この場合には P''^t -a の正となる要素のいずれも固定することが出来ますが、
全体を一度に固定することは出来ません。x:=P''^t -a を先に計算し、
||a|| を適切に調整して再度計算し直すことによるループによって
通常はきちんと計算出来ますが、誤差に関する変動がとてもシビアです。
これは主にベクトルの射影の要素の符号の誤差によるものです。

また、このアルゴリズムは最終段の切片のスケーリングにより原点周辺で
精度を著しく毀損します。

# 使い方
    // 必要であれば namespace ブロックで括ってください。ただし、インクルードガードが有害な場合があります。
    #include "konbu_init.h"
    ...
    // mpfr の場合, num_t::set_default_proc(BITS); が必要です
    // QD   の場合,   unsigned int old_cw; fpu_fix_start(&old_cw); が必要です
    ...
    Mat A(/* some row */, /* some col */);
    Vec b(A.rows());
    ...
    const auto error(A * Linner<num_t>().inner(A, b) - b);
    ...
    // QD   の場合,    fpu_fix_end(&old_cw); が必要です

# 少しの情報
もし、mn コア以上のコア数で実行している場合には、O(n\*lg m) 時間で内点を取ってこれるように書き換えることができます。  
また、もし一つ一つ固定せずに、一度に次元の数だけ固定できれば、最終的に O(lg(n)\*lg(mn)) 時間で実行できるようになる可能性がありますが、望み薄です。これは、不要な超平面が含まれる可能性があるためです。

# その他のダウンロードサイト
* https://konbu.azurewebsites.net/ (Sample Site)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/
* https://www.sourceforge.net/projects/convcheck/

# アーカイブ
このレポジトリはアーカイブされました。バグレポート以外では変更の予定はありません。

