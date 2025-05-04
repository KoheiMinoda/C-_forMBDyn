// 標準ライブラリ

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// <vector> を使ったデータ管理クラス
class DataList {
    private:
        // DataList の中で使う変数の定義
        std::vector<double> data;

    public:
        void add(double value) {
            data.push_back(value);
        }

        const std::vector<double>& getData() const {
            return data;
        }

        void print() const {
            std::cout << "DataList: "; // std::endl; が無い → 改行もしない
            // data の中身の一つ一つを "d" とし，中のコードを実行する
            for (double d : data) {
                std::cout << d << " "; // std::endl; が無い → 改行もしない
            }
            std::cout << std::endl;
        }
};

// <cmath> を使った数学ツールクラス
class MathTools {
    public:
        static double getMagnitude(double x) {
            return std::abs(x);
        }

        static double getSqrt(double x) {
            return std::sqrt(x);
        }

        static double getSin(double rad) {
            return std::sin(rad);
        }
};

// <algorithm> を使ったデータ処理クラス
class DataProcessor {
    public:
        static double getMax(const std::vector<double>& vec) {
            // max_element は　algorithm　ライブラリに搭載されている機能
            return *std::max_element(vec.begin(), vec.end());
        }

        static void sortAscending(std::vector<double>& vec) {
            // sort は　algorithm　ライブラリに搭載されている機能
            std::sort(vec.begin(), vec.end());
        }
};


// main()：各クラスを使って動作確認
int main() {
    // vector機能
    DataList list;
    list.add(3.0);
    list.add(-1.5);
    list.add(4.2);
    list.print();

    // cmath機能
    double x = -2.5;
    // getMagnitude, getSqrt, getSin はそれぞれ MathTools ライブラリに搭載されている機能
    std::cout << "絶対値: " << MathTools::getMagnitude(x) << std::endl;
    std::cout << "平方根: " << MathTools::getSqrt(9.0) << std::endl;
    std::cout << "sin(π/2): " << MathTools::getSin(M_PI / 2) << std::endl;

    // algorithm機能
    std::vector<double> data = list.getData();  // DataList list; をここでも使って小数格納ベクトル "data" を作成する
    double maxVal = DataProcessor::getMax(data);
    DataProcessor::sortAscending(data);

    std::cout << "最大値: " << maxVal << std::endl;
    std::cout << "ソート後: ";
    for (double d : data) std::cout << d << " ";
    std::cout << std::endl;

    return 0;
}
