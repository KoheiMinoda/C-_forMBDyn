// 繰り返し

#include <iostream> 
#include <string>

void simpleMotion(double v, double g, double dt) {
    // 初期値; 範囲; 更新
    for (int t=0; t<5; ++t) {
        v += g * dt;
        std::cout << t+1 << "秒後の速度は" << v << "m/s" << std::endl;
    }
}

void simpleMotion2(double u, double a, double delta_t) {
    std::cout << "-------------" << std::endl;
    // 初期値
    int t=0;
    // 範囲
    while (t<=10) {
        u += a * delta_t;
        std::cout << t+1 << "秒後の速度は" << u << "m/s" << std::endl;
        //更新
        t++;
    }   
}

void evenSkipCount(int max) {
    std::cout << "-------------" << std::endl;
    int sum = 0;
    for (int i=0; i<max; i++) {
        // 30 以上はやらない
        if (i == 30) break;
        // 偶数はスキップ
        if (i % 2 == 0) continue;
        // 奇数だけ積み上げていく
        sum += i;
    }
    std::cout << "奇数の合計は: " << sum << std::endl;
}

int main() {

    double v=0.0; //[m/s]
    double g=9.8; // [m/s^2]
    double dt=1; // [t]
    simpleMotion(v, g, dt);

    double u=20.0; //[m/s]
    double a=-2.75; // [m/s^2]
    double delta_t=0.2; // [t]
    simpleMotion2(u, a, delta_t);

    int max=20;
    evenSkipCount(max);

    return 0;
}