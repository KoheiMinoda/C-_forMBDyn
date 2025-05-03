// 条件分岐

#include <iostream> 
#include <string> 

double force = 0.0;
bool isFixed = true;

// if 条件分岐
void checkConstraint(bool isFixed) {
    if (isFixed) {
        std::cout << "固定" << std::endl;
    } else {
        std::cout << "移動可能" << std::endl;
    }
}

// switch
void judgeMode(int mode) {
    switch (mode) {
        case 1:
            std::cout << "mode 1" << std::endl;
            break;
        case 2:
            std::cout << "mode 2" << std::endl;
            break;
        case 3:
            std::cout << "mode 3" << std::endl; 
    }
}

// main 関数の中に詰め込むよりも，外で関数を定義して main 関数は実行に専念させることで可読性が良くなる
int main() {
    // nodeFixed が false の場合の checkConstraint 関数の実行
    bool nodeFixed = false;
    checkConstraint(nodeFixed);

    // mode が 2 の場合の judgeMode 関数の実行
    int mode = 2;
    judgeMode(mode);

    return 0;
}