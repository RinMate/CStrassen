StrassenのアルゴリズムをC言語で実装しました。
サンプルコードでは一つの行列を読み込みそれをk乗するプログラムを実装しています。

## 入力
`input.txt` にあるように  
1行目：行と列の数を整数値で与えます n m  
2~n+1行目：行列の中身を数値で与えます  
n+2行目：累乗する回数を与えます k

## 出力
入力と同じような形式で出力します
```
$ ./a.out < input.txt
5 5
1329.3450000000 445.2500000000 1434.3000000000 1297.9000000000 1491.8500000000
1662.7200000000 459.7210000000 1388.1000000000 1568.2000000000 1665.0700000000
934.2000000000 571.1300000000 1355.0450000000 1001.7800000000 1138.2200000000
994.2600000000 251.8500000000 869.4000000000 930.6450000000 863.7000000000
2908.6900000000 943.4200000000 3058.3000000000 3209.0500000000 2979.3000000000
```