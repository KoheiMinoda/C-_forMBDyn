// 配列

#include <iostream>

/*

基本的な使い方

#include <iostream>

int main() {
    double velocities[5] = {0.0, 1.0, 2.0, 3.0, 4.0};

    for (int i = 0; i < 5; ++i) {
        std::cout << i << " 番目の速度は " << velocities[i] << " m/s" << std::endl;
    }

    return 0;
}

*/

// 位置と速度の情報を持たせる
int main() {
    const int N = 3; // 3 つのノードについて考える
    double position[N] = {0.0, 1.0, 2.0}; // 三次元の情報を持たせるためにはここで工夫が必要
    double velocity[N] = {0.0, 0.5, 1.0}; // 三次元の情報を持たせるためにはここで工夫が必要
    double dt = 1.0;

    for (int i = 0; i < N; ++i) {
        position[i] += velocity[i] * dt;
        std::cout << "ノード" << i << ": 位置 = " << position[i] << std::endl;
    }

    return 0;
}

