// 構造体とクラス

#include <iostream>
#include <string>

// struct : 構造体　：　データのまとまり
// 変数や関数を保存しているが，これ単体で実行されることはなく，主に関数にデータを渡す役割となる
struct Node {
    double x;  // 位置
    double v;  // 速度
    double a;  // 加速度

    void update(double dt) {
        v += a * dt;
        x += v * dt;
    }

    void print(int index) {
        std::cout << "ノード" << index << " → 位置: " << x << " m, 速度: " << v << " m/s" << std::endl;
    }
};

// class : 構造体 + カプセル化・アクセス制御
class Simulation {
    // private: アクセス制御 → 堅牢なコードに繋がる
    // 例えば，下の private で定義している nodes と dt を main() 関数の中で更新できない
    // Simulation クラスの中でだけ使える
    private:
        // 上で定義した struct Node の情報を使う
        Node nodes[3]; // ３つのノードについて考える
        double dt;

    public:
        // public: の中の関数には，main() 関数の中で引数を渡したりできる
        // コンストラクタ → オブジェクトが生成されるときに最初に1回だけ自動的に呼ばれる関数（詳しくは後述）
        Simulation(double delta_t) {
            dt = delta_t;

            // それぞれのノードの初期化：位置, 初速度, 加速度
            nodes[0] = {0.0, 1.0, 0.5};    // 加速するノード
            nodes[1] = {2.0, 0.0, -1.0};   // 減速するノード
            nodes[2] = {5.0, -0.5, 0.0};   // 等速直線運動
        }

        // メソッド → クラスや構造体の中に定義された関数 : オブジェクトに対して操作を行う
        // 関数 → クラス外にある普通の処理ブロック
        void runStep() {
            for (int i = 0; i < 3; ++i) {
                nodes[i].update(dt);
                nodes[i].print(i);
            }
        }
};


int main() {
    double dt = 1.0;  // Simulation 関数を実行するために必要な変数を定義

    // Simulation クラスの"オブジェクト（インスタンス）" sim を作成し，dt という値で初期化する
    Simulation sim(dt);

    std::cout << "===== 1ステップ後の状態 =====" << std::endl;
    // 公開されたメソッドを実行
    sim.runStep();

    return 0;
}
