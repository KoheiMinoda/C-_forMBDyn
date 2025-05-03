// 変数とデータ型

#include <iostream>
#include <string>  // 文字列を使うとき必要：01 では C 言語由来の文字列をそのまま出力するので不要

int main() {
    double mass = 10.5;
    double acceleration = 9.8;
    double force = mass * acceleration;

    std::cout << "質量: " << mass << " kg" << std::endl;
    std::cout << "加速度: " << acceleration << " m/s^2" << std::endl;
    std::cout << "力: " << force << " N" << std::endl;

    return 0;
}

/*

MBDyn との関係

double 型：小数：剛体の位置・速度・加速度・力などで頻出
bool型：TRUE, FALSE：拘束条件や「固定されているか」などの判定に使用
string 型：文字列：部品ID、ノード名、構造体のキーなど

*/