// 関数

#include <iostream> 
#include <string>

// void ：　戻り値を有さない関数の宣言
void greet(std::string name) {
    std::cout << "こんにちは" << name << "さん！" << std::endl;
}

// double : 小数の戻り値を返す関数の宣言
double calcForce(double mass, double acc) {
    return mass * acc;
}

// int : 整数の戻り値を返す関数の宣言
// main() は事前に C++ の仕様で特別に定義されている"エントリーポイント"（ここでプログラムの実行が始まる）
int main() {
    greet("kohei");

    double F = calcForce(2.0, 9.8);
    std::cout << "力: " << F << " N" << std::endl;

    // これまでなんとなくつけていた main() 関数は， return 0; で整数 0 を返すということで，OS に正常終了を伝えている
    // int main() { return 0; } は C++ というプログラミング言語において特別な意味を持っている
    return 0;
}

/*

MBDyn内で使う User Defined Element (UDE) では main() を書かずに，指定された名前の関数（例：UserElemInit()）がエントリーポイントになる

ただし，main() を通じて理解した概念（戻り値・引数・構文など）は UDE 設計に応用可能

*/