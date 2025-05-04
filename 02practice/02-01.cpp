// クラスの継承

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

// 基底クラス
// 基底クラス単体ではインスタンス化できず，あくまで土台として使う
class Element {
    // protected: 派生クラスからはアクセス可能
    protected:
        std::string name;
    public:
        // コンストラクタで名前を初期化
        Element(const std::string& n) : name(n) {}
        // virtual は仮想関数
        // Element の派生クラスでは，必ず update と print を定義しなければならない
        virtual void update(double dt) = 0;
        virtual void print() const = 0; // const はメンバ変数を書き換えないという宣言（安全な出力に繋がる）
        virtual ~Element() {} 
};

class RigidBody : public Element {
    private:
        // 位置 x と速度 v を持つ剛体
        double x, v;
    public:
        // コンストラクタで名前，初期位置，初速度を設定
        RigidBody(std::string n, double x0, double v0)
            : Element(n), x(x0), v(v0) {}

        void update(double dt) override {
            x += v * dt;
        }

        void print() const override {
            std::cout << name << " [RigidBody] 位置: " << x << ", 速度: " << v << std::endl;
        }
};

class Spring : public Element {
    private:
        double k, x;
    public:
        Spring(std::string n, double stiffness, double x0)
            : Element(n), k(stiffness), x(x0) {}

        void update(double dt) override {
            double force = -k * x;
            x += force * dt;
        }

        void print() const override {
            std::cout << name << " [Spring]     変位: " << x << ", 力: " << -k * x << std::endl;
        }
};

class Damper : public Element {
    private:
        double c, v;
    public:
        Damper(std::string n, double damping, double v0)
            : Element(n), c(damping), v(v0) {}

        void update(double dt) override {
            double force = -c * v;
            v += force * dt;
        }

        void print() const override {
            std::cout << name << " [Damper]     速度: " << v << ", 力: " << -c * v << std::endl;
        }
};

int main() {
    std::vector<Element*> elements;
    elements.push_back(new RigidBody("Body1", 0.0, 1.0));
    elements.push_back(new Spring("Spring1", 10.0, 0.5));
    elements.push_back(new Damper("Damper1", 5.0, 2.0));

    const int steps = 10;
    const double dt = 0.1;

    for (int t = 0; t <= steps; ++t) {
        std::cout << "\n--- 時刻: " << std::fixed << std::setprecision(1) << t * dt << " s ---" << std::endl;
        for (auto* e : elements) {
            e->update(dt);
            e->print();
        }
    }

    // メモリ解放
    for (auto* e : elements) delete e;

    return 0;
}
