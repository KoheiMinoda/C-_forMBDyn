#include "mbconfig.h"           // This goes first in every *.c,*.cc file

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <limits>

#include "dataman.h"
#include "userelem.h"
#include "module-catenary_lm.h" // 自前のヘッダ

// 追加のヘッダー
// dataman.h や userelem.h から間接的にインクルードされるならいらない可能性
#include "elem.h"          // Elem 基底クラス (userelem.hがインクルードしている可能性が高い)
#include "strnode.h"       // StructNode (ランプドマス点の表現に必須)
#include "drive.h"         // DriveOwner (FSFのため)
#include "node.h"          // Node (pGetNode の戻り値の基底クラス)
#include "gravity.h"

// 追加の標準ライブラリ
#include <vector>       // std::vector (ノードやセグメントのリスト管理用)
#include <numeric>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159
#endif

// クラス宣言
class ModuleCatenaryLM : virtual public Elem, public UserDefinedElem
{
public:
    // コンストラクタ [uLabel : 要素のラベル, *pD0 : DOF(Degree of Freedom) 所有者へのポインタ, *pDM : データマネージャーへのポインタ, HP : MBDyn パーサへの"参照"] 
    ModuleCatenaryLM(unsigned uLabel, const DofOwner *pD0, DataManager* pDM, MBDynParser& HP);
    // デストラクタ [オブジェクトが破棄されるたび呼びだされる]
    // メモリリークを防ぐために重要
    virtual ~ModuleCatenaryLM(void);

    // ========= 仮想関数群 ここから ==============
    // クラス宣言の際に Elem と UserDefinedElem を継承している
    // 以下の "virtual" 宣言は (Elem と UserDefinedElem) で宣言されている仮想関数を再定義することを意味する
    // "const" のキーワードがついているメンバ関数は，その関数内でオブジェクトのメンバ変数を変更しないという宣言

    // シミュレーション中に特定の要素の状態（FP 点の位置や張力など）を指定された形式で出力ファイルに書き出すために MBDyn から呼び出される
    virtual void Output(OutputHandler& OH) const;
    // ヤコビアン行列，残差行列のサイズを決める：ランプドマスモデルでは，関連するノードの自由度に依存
    virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

    // 要素の運動方程式を現在の状態変数とその時間微分で偏微分下ヤコビアン行列の"成分"を計算し，MBDyn の全体ヤコビアン行列の対応する部分に格納する
    VariableSubMatrixHandler&
    AssJac(VariableSubMatrixHandler& WorkMat, 
        doublereal dCoef, 
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr
    );

    // 要素の運動方程式の外力を計算する
    SubVectorHandler&
    AssRes(SubVectorHandler& WorkVec, 
        doublereal dCoef, 
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr
    );

    // MBDyn ソルバーがノードの変位や速度を更新した後，要素の内部状態をそれに応じて更新するために呼び出す
    virtual void SetValue(DataManager* pDM,
        VectorHandler& X,
        VectorHandler& XP,
        SimulationEntity::Hints* pHints
    ) override; // override は C++11 以降で基底クラスに仮想関数がなければコンパイルエラーを発生させるという機能：無くても良い

    // この要素が直接族または管理しているノードの数を返す
    virtual unsigned int iGetNumConnectedNodes(void) const;
    // シミュレーション開始時に要素の初期状態や初期ノードの初期値を設定する
    virtual void SetInitialValue(VectorHandler& X, VectorHandler& XP);
    // 中断した後に再開する場合に内部情報を保存する
    virtual std::ostream& Restart(std::ostream& out) const;

    virtual unsigned int iGetInitialNumDof(void) const;
    virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    virtual VariableSubMatrixHandler&
    InitialAssJac(VariableSubMatrixHandler& WorkMat, const VectorHandler& XCurr);
    virtual SubVectorHandler&
    InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

    virtual const Node* pGetNode(unsigned int i) const;
    virtual Node* pGetNode(unsigned int i);
    // ========= 仮想関数群 ここまで ==============

private:
    double APx_orig; // アンカー点の x 座標
    double APy_orig; // アンカー点の y 座標
    double APz_orig; // アンカー点の z 座標
    double L_orig; // ライン全長
    double w_orig; // チェーンの単重
    double xacc_orig; // rtsafe 関数で使用するパラメータ：収束判定

    DriveOwner FSF_orig; // 元の Force Scale Factor：ランプアップに使う

    // ===== 静的メンバ関数 =====
    // "static" 特定のオブジェクトインスタンスに属さず，クラス自体に関連付けられる
    static double myasinh_local(double val);
    static double myacosh_local(double val);
    static double myatanh_local(double val);
    static void funcd_catenary_local(double x_param, double xacc, double &f_val, double &df_val, double d_geom, double l_geom, double &p0_angle);
    static double rtsafe_catenary_local(double x1_bounds, double x2_bounds, double xacc_tol, double d_geom, double l_geom, double &p0_angle_out);
    // ========================

    unsigned int Seg_param;

    // ノード配列 (StructNode* の動的配列)
    // N_nodes_param[0] はフェアリーダーノード
    // N_nodes_param[1] ... N_nodes_param[Seg_param-1] は内部質量点ノード
    // アンカーは固定座標(APx_orig, APy_orig, APz_orig)として別途扱う
    // このベクトルサイズは Seg_param になる（0 ~ Seg_param - 1）
    std::vector<StructDispNode*> N_nodes_param;

    // 各セグメントの特性
    struct SegmentProperty {
        doublereal L0_seg;  // 「セグメント」の初期自然長
        doublereal M_seg;   // 「セグメント」の質量
        doublereal EA_seg;  // 「セグメント」の軸剛性 (EA)
        doublereal CA_seg;  // 「セグメント」の軸方向減衰係数 (オプション)

        SegmentProperty() : L0_seg(0.0), M_seg(0.0), EA_seg(0.0), CA_seg(0.0) {}
    };
    std::vector<SegmentProperty> P_param; // サイズは Seg_param (セグメント数)

    doublereal g_gravity_param;
};

// コンストラクタ：パラメータの読み込みと初期設定
ModuleCatenaryLM::ModuleCatenaryLM(
    unsigned uLabel,
    const DofOwner *pD0,
    DataManager* pDM,
    MBDynParser& HP
)
    // 初期化子リスト
    : Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pD0),
        // メンバ変数の初期化
        APx_orig(0.0), APy_orig(0.0), APz_orig(0.0),
        L_orig(0.0), w_orig(0.0), xacc_orig(1e-6),
        Seg_param(20),
        g_gravity_param(9.80665)
    {
        // 入力ファイルの現在位置が "help" なら使用方法のメッセージを表示
        if (HP.IsKeyWord("help")) {
            silent_cout(
                "\n"
                "Module:     ModuleCatenaryLM (Lumped Mass Catenary Mooring Line - Segments Fixed to 20)\n"
                "Usage:      catenary_lm, label,\n"
                "                fairlead_node_label,\n"
                "                total_length, unit_weight, rtsafe_accuracy,\n"
                "                APx, APy, APz, (Anchor fixed coordinates)\n"
                "                EA (axial_stiffness),\n"
                "              [ CA (axial_damping), ]\n"
                "              [ gravity, g, ]\n"
                "              [ force scale factor, (DriveCaller), ]\n"
                "              [ output, (FLAG) ] ;\n"
                "\n"
                << std::endl
            );
            if (!HP.IsArg()) {
                throw NoErr(MBDYN_EXCEPT_ARGS);
            }
        }

        // ===== パラメータ読み込み開始 =====
        // ここでの読み込みは .mbd ファイルの elements block を上から読み込む（.usr ファイルに外だししているかも）

        // 全長 L_orig
        L_orig = HP.GetReal();
        if (L_orig <= 0.0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        // 単位重量 w_orig
        w_orig = HP.GetReal();
        if (w_orig <= 0.0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        // rtsafeの計算精度 xacc_orig
        xacc_orig = HP.GetReal();
        if (xacc_orig <= 0.0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        // アンカー固定座標 APx_orig, APy_orig, APz_orig
        APx_orig = HP.GetReal();
        APy_orig = HP.GetReal();
        APz_orig = HP.GetReal();

        // ====== EA, CA, g : 既存の .usr ファイルに加える項目？==========
        // 軸剛性 EA (ランプドマス法用)
        doublereal EA_val;
        if (HP.IsKeyWord("EA")) { // キーワードがある方がより頑健
            EA_val = HP.GetReal();
        } else {
            EA_val = HP.GetReal(); // キーワードなしで順番に読む場合
        }
        if (EA_val <= 0.0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        // 軸方向減衰係数 CA (オプション)
        doublereal CA_val = 0.0; // デフォルトは減衰なし
        if (HP.IsKeyWord("CA")) {
            CA_val = HP.GetReal();
            if (CA_val < 0.0) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        // 重力加速度 g (オプション)
        if (HP.IsKeyWord("gravity")) {
            g_gravity_param = HP.GetReal();
            if (g_gravity_param < 0.0) { // 正の値
            }
        }

        // FSF の読み込み（元のコードを維持）
        if (HP.IsKeyWord("Force" "scale" "factor")) {
			FSF_orig.Set(HP.GetDriveCaller());
		} else {
			FSF_orig.Set(new OneDriveCaller);
		}

        // 出力フラグの設定（元のコードを維持）
        SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

        // ノードの確保と初期化
        // N_nodes_param にはフェアリーダーと Seg_param-1 個の内部ノードを格納
        // 合計 Seg_param 個の StructNode* を持つ
        N_nodes_param.clear();
        N_nodes_param.resize(Seg_param);

        // ====== ここも既存の .usr ファイルに書き足す (201~220 などユニークな数字を 20 コ) ======
        unsigned int fairlead_node_label = HP.GetInt();

        // 各ノードのラベルの受け取り
        // フェアリーダーノードの処理
        {
            Node* rawNode = pDM->pFindNode(Node::STRUCTURAL, fairlead_node_label);
            if (rawNode == nullptr) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
            StructDispNode* dispNode = dynamic_cast<StructDispNode*>(rawNode);
            if (dispNode == nullptr) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
            N_nodes_param[0] = dispNode;
        }

        // 内部質量点ノードの処理（1 … Seg_param-1）もすべて入力ファイルでラベルを受け取り
        for (unsigned int i = 1; i < Seg_param; ++i) {
            unsigned int node_label = HP.GetInt();

            Node* rawNode = pDM->pFindNode(Node::STRUCTURAL, node_label);
            if (rawNode == nullptr) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
            StructDispNode* dispNode = dynamic_cast<StructDispNode*>(rawNode);
            if (dispNode == nullptr) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
            N_nodes_param[i] = dispNode;
        }

        Vec3 zeroPos(0.0, 0.0, 0.0);
        Vec3 zeroVel(0.0, 0.0, 0.0);

        P_param.resize(Seg_param); // Seg_param 個のセグメント
        doublereal L0_s = L_orig / static_cast<doublereal>(Seg_param); // 各セグメントの自然長：単純に全長をセグメント数で割っている
        doublereal mass_per_segment = (w_orig / g_gravity_param) * L0_s; // 各セグメントの質量：各セグメントの「質量」

        // P_param ベクトルの各要素にそれぞれを格納
        for (unsigned int i = 0; i < Seg_param; ++i) {
            P_param[i].L0_seg = L0_s;
            P_param[i].M_seg = mass_per_segment;
            P_param[i].EA_seg = EA_val;
            P_param[i].CA_seg = CA_val;
        }

        // ログ出力：
        pDM->GetLogFile() << "ModuleCatenaryLM (" << GetLabel() << ") initialized:" << std::endl;
        pDM->GetLogFile() << "  Fairlead Node Label: " << fairlead_node_label << std::endl;
        pDM->GetLogFile() << "  Anchor Fixed At: (" << APx_orig << ", " << APy_orig << ", " << APz_orig << ")" << std::endl;
        pDM->GetLogFile() << "  Original Line Length: " << L_orig << ", Unit Weight: " << w_orig << std::endl;
        pDM->GetLogFile() << "  Segments (Fixed): " << Seg_param << std::endl;
        pDM->GetLogFile() << "  EA: " << EA_val << ", CA: " << CA_val << ", Gravity: " << g_gravity_param << std::endl;
        pDM->GetLogFile() << "  RTSAFE Accuracy: " << xacc_orig << std::endl;
    }

ModuleCatenaryLM::~ModuleCatenaryLM(void){}

double ModuleCatenaryLM::myasinh_local(double val) { return std::log(val + std::sqrt(val * val + 1.)); }
double ModuleCatenaryLM::myacosh_local(double val) { return std::log(val + std::sqrt(val + 1.) * std::sqrt(val - 1.));}
double ModuleCatenaryLM::myatanh_local(double val) { return 0.5 * std::log((1. + val) / (1. - val)); }

// rtsafe_catenary_local 関数関数から呼び出される：f(x) とその導関数を求める
void ModuleCatenaryLM::funcd_catenary_local(double x_param, double xacc, double &f_val, double &df_val, double d_geom, double l_geom, double &p0_angle) {
    int i_internal_iter, max_internal_iter;
    double f1_internal, df1_internal;
    max_internal_iter = 1000;

    // 水平張力が無い状態
    if(x_param == 0.0) {
        f_val = -d_geom;
        df_val = 0.0;
        p0_angle = 0.0;

    // 水平張力がある状態
    } else if(x_param > 0.0) {
        
        // 特殊なケース：元のコードをそのままにしつつ，保護をすこししただけ
        if(l_geom <= 0.0){
            double X_1_internal;
            X_1_internal = 1.0/x_param+1.0;
            if (X_1_internal < 1.0) X_1_internal = 1.0; // acoshの引数保護
                double term_sqrt1 = 1.0+2.0*x_param;
            if (term_sqrt1 < 0.0) term_sqrt1 = 0.0;
                double term_sqrt2 = X_1_internal*X_1_internal-1.0;
            if (term_sqrt2 < 0.0) term_sqrt2 = 0.0;

            f_val=x_param*myacosh_local(X_1_internal)-std::sqrt(term_sqrt1)+1.0-d_geom;
            if (std::fabs(x_param * std::sqrt(term_sqrt2)) < 1e-12 || std::fabs(std::sqrt(term_sqrt1)) < 1e-12 ) { // ゼロ除算保護
                df_val = myacosh_local(X_1_internal); // 近似
            } else {
                df_val=myacosh_local(X_1_internal)-1.0/std::sqrt(term_sqrt1)-1.0/(x_param*std::sqrt(term_sqrt2));
            }
            p0_angle=0.0;

        // 一般的なケース
        } else {

            // 海底との接触がある可能性がある場合：元のコードをそのままにしつつ，保護をすこししただけ
            if(x_param > (l_geom*l_geom - 1.0) / 2.0) {
                p0_angle=0.0;
                for(i_internal_iter = 0; i_internal_iter < max_internal_iter; i_internal_iter++) {
                    double cos_p0 = std::cos(p0_angle);
                    if (std::fabs(cos_p0) < 1e-9) { df1_internal = 1.0; break; } // 保護
                    double func1_internal = 1.0/x_param + 1.0/cos_p0;
                    double term_in_sqrt_f1 = func1_internal*func1_internal - 1.0;
                    if (term_in_sqrt_f1 < 0.0) term_in_sqrt_f1 = 0.0;

                    f1_internal = x_param*(std::sqrt(term_in_sqrt_f1) - std::tan(p0_angle)) - l_geom;

                    if (std::fabs(cos_p0) < 1e-9 || term_in_sqrt_f1 < 1e-12 ) { df1_internal = 1.0; break; } // 保護
                    df1_internal = x_param * (func1_internal * std::tan(p0_angle) / (cos_p0 * std::sqrt(term_in_sqrt_f1)) - (std::tan(p0_angle)*std::tan(p0_angle)) - 1.0);

                    if (std::fabs(df1_internal) < 1e-9) { break; }
                    p0_angle = p0_angle-f1_internal/df1_internal;

                    cos_p0 = std::cos(p0_angle);
                    if (std::fabs(cos_p0) < 1e-9) { break; }
                    func1_internal = 1.0/x_param + 1.0/cos_p0;
                    term_in_sqrt_f1 = func1_internal*func1_internal - 1.0;
                    if (term_in_sqrt_f1 < 0.0) term_in_sqrt_f1 = 0.0;
                    f1_internal = x_param*(std::sqrt(term_in_sqrt_f1) - std::tan(p0_angle)) - l_geom;

                    if(std::fabs(f1_internal) < xacc) { break; }
                }
                if(i_internal_iter == max_internal_iter && std::fabs(f1_internal) > xacc) {
                }

                double X_2_internal = l_geom/x_param + std::tan(p0_angle);
                double X_3_internal = std::tan(p0_angle);
                f_val = x_param*(myasinh_local(X_2_internal) - myasinh_local(X_3_internal)) - l_geom + 1.0 - d_geom;
                
                double term_in_sqrt_df = X_2_internal*X_2_internal + 1.0;
                // if (term_in_sqrt_df < 0.0) term_in_sqrt_df = 0.0; // asinhの引数は実数なので常に正
                if (std::fabs(x_param * std::sqrt(term_in_sqrt_df)) < 1e-12 ) { df_val = 1.0; } // 保護
                else {
                    df_val=myasinh_local(X_2_internal) - myasinh_local(X_3_internal) - l_geom/(x_param*std::sqrt(term_in_sqrt_df));
                }

            // 海底との接触がある可能性が無い場合：単純なカテナリー
            } else {
                double X_5_internal = 1.0/x_param+1.0;
                if (X_5_internal < 1.0) X_5_internal = 1.0; // acoshの引数保護
                double term_sqrt1 = 1.0+2.0*x_param;
                if (term_sqrt1 < 0.0) term_sqrt1 = 0.0;
                double term_sqrt2 = X_5_internal*X_5_internal-1.0;
                if (term_sqrt2 < 0.0) term_sqrt2 = 0.0;

                f_val = x_param*myacosh_local(X_5_internal) - std::sqrt(term_sqrt1) + 1.0 - d_geom;
                if (std::fabs(x_param * std::sqrt(term_sqrt2)) < 1e-12 || std::fabs(std::sqrt(term_sqrt1)) < 1e-12) { // 保護
                     df_val = myacosh_local(X_5_internal); // 近似
                } else {
                     df_val = myacosh_local(X_5_internal) - 1.0/std::sqrt(term_sqrt1) - 1.0/(x_param*std::sqrt(term_sqrt2));
                }
                p0_angle = 0.0;
            }
        }
    } else {
        throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
    }
}

// x1_bounds, x2_bounds は根が含まれると期待される初期区間の下限と上限
double ModuleCatenaryLM::rtsafe_catenary_local(double x1_bounds, double x2_bounds, double xacc_tol, double d_geom, double l_geom, double &p0_angle_out) {
    const int MAXIT_internal = 1000;
    int j_internal;
    double fh_internal,fl_internal,xh_internal,xl_internal;
    double dx_internal,dxold_internal,f_internal,temp_internal,rts_internal;
    double p1_internal, p2_internal;
    double df_not_used, df_internal; // funcd_catenary_local がdfを要求するため

    // x1_bounds, x2_bounds における f(x) とその導関数
    funcd_catenary_local(x1_bounds, xacc_tol, fl_internal, df_not_used, d_geom, l_geom, p1_internal);
    funcd_catenary_local(x2_bounds, xacc_tol, fh_internal, df_not_used, d_geom, l_geom, p2_internal);

    // 二つの関数値が同符号ならエラー
    if((fl_internal > 0.0 && fh_internal > 0.0) || (fl_internal < 0.0 && fh_internal < 0.0)) {
        throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
    }
    // どちらかの関数値が 0 であればそれが答え
    if(fl_internal == 0.0) {
        p0_angle_out = p1_internal;
        return x1_bounds;
    }
    if(fh_internal == 0.0) {
        p0_angle_out = p2_internal;
        return x2_bounds;
    }

    // fl_internal < 0 となるように xl_internal と xh_internal を設定する
    if(fl_internal < 0.0) {
        xl_internal = x1_bounds;
        xh_internal = x2_bounds;
    } else {
        xh_internal = x1_bounds;
        xl_internal = x2_bounds;
    }

    // 中点を初期の推定値とし，反復処理のセット
    rts_internal = 0.5*(x1_bounds + x2_bounds);
    dxold_internal = std::fabs(x2_bounds - x1_bounds);
    dx_internal = dxold_internal;
    funcd_catenary_local(rts_internal, xacc_tol, f_internal, df_internal, d_geom, l_geom, p0_angle_out);


    for(j_internal = 0; j_internal < MAXIT_internal; j_internal++) {
        if((((rts_internal - xh_internal)*df_internal - f_internal)*((rts_internal - xl_internal)*df_internal - f_internal) > 0.0)
           || (std::fabs(2.0*f_internal) > std::fabs(dxold_internal*df_internal))) {
            dxold_internal = dx_internal;
            dx_internal = 0.5*(xh_internal - xl_internal);
            rts_internal = xl_internal + dx_internal;
            if(xl_internal == rts_internal) { return rts_internal; }
        } else {
            dxold_internal = dx_internal;
            if (std::fabs(df_internal) < 1e-12) { // ゼロ除算を避ける
                dx_internal = 0.5*(xh_internal - xl_internal); // 二分法にフォールバック
                rts_internal = xl_internal + dx_internal;
                if(xl_internal == rts_internal){ return rts_internal; }
            } else {
                dx_internal = f_internal/df_internal;
            }
            temp_internal = rts_internal;
            rts_internal -= dx_internal;
            if(temp_internal == rts_internal) {return rts_internal;}
        }

        if(std::fabs(dx_internal) < xacc_tol) { return rts_internal; }

        funcd_catenary_local(rts_internal, xacc_tol, f_internal, df_internal, d_geom, l_geom, p0_angle_out);

        if(f_internal < 0.0){
            xl_internal = rts_internal;
        } else {
            xh_internal = rts_internal;
        }
    }

    throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
}

static inline double local_robust_acosh_arg(double val) {
    const double min_arg = 1.0;
    return (val < min_arg) ? min_arg : val;
}

// 要素に関連する自由度の初期値を設定するために MBDyn から呼び出す
// N_nodes_params に格納されているノードの初期位置や初期速度を計算
void ModuleCatenaryLM::SetInitialValue(VectorHandler& /* X */, VectorHandler& /* XP */) {
    // フェアリーダー点座標とアンカー点座標からラインの初期形状を計算し，各内部ノードの初期座標を決定する，かも？
    if (N_nodes_param.empty() || N_nodes_param[0] == 0) {
        silent_cerr("ModuleCatenaryLM::SetInitialValue: ERROR: Fairlead node (N_nodes_param[0]) is not set or null." << std::endl);
        throw ErrGeneric("ModuleCatenaryLM::SetInitialValue: Fairlead node not properly initialized.");
    }
    if (w_orig < std::numeric_limits<double>::epsilon()) {
        silent_cerr("ModuleCatenaryLM::SetInitialValue: WARNING: Unit weight w_orig is zero or very small. Catenary calculations may be unstable or lead to a straight (taut) line if L_orig matches direct distance." << std::endl);
    }

    const Vec3 A_pos_global(APx_orig, APy_orig, APz_orig);
    const Vec3 F_pos_global = N_nodes_param[0]->GetXCurr().pos;

    const Vec3 AF_vec_global = F_pos_global - A_pos_global;
    const double L_APFP_sq = AF_vec_global.dGet(1)*AF_vec_global.dGet(1) + AF_vec_global.dGet(2)*AF_vec_global.dGet(2);
    const double L_APFP = std::sqrt(L_APFP_sq);
    const double h_FP_rel_AP = AF_vec_global.dGet(3);
    const h_abs = std::fabs(h_FP_rel_AP);

    const double epsilon = 1.e-9;

    const double L_straight_sq = L_APFP_sq + h_FP_rel_AP * h_FP_rel_AP;
    const double L_straight = std::sqrt(L_straight_sq);

    if (L_orig < L_straight - epsilon) {
        silent_cerr("ModuleCatenaryLM::SetInitialValue: ERROR: Line length L_orig (" << L_orig
                  << ") is less than the straight distance between anchor and fairlead (" << L_straight << ")." << std::endl);
        for (unsigned int i = 0; i < Seg_param; ++i) {
            if (N_nodes_param[i]) {
                N_nodes_param[i]->XCurr.pos = F_pos_global;
                N_nodes_param[i]->VCurr.v = Vec3(0., 0., 0.);
            }
        }
        throw ErrGeneric("ModuleCatenaryLM::SetInitialValue: Line too short.");
    }

    if (std::fabs(L_orig - L_straight) < epsilon || w_orig < epsilon ) {
        if (w_orig < epsilon && L_orig > L_straight + epsilon) {
             silent_cerr("ModuleCatenaryLM::SetInitialValue: WARNING: Line is weightless and longer than straight distance. Assuming taut configuration for initial setup." << std::endl);
        }

        for (unsigned int k_node_idx = 0; k_node_idx < Seg_param; ++k_node_idx) {
            double s_arc_ratio = static_cast<double>(Seg_param - k_node_idx) / static_cast<double>(Seg_param);
            double s_arc_from_anchor = L_orig * s_arc_ratio;
            
            Vec3 node_pos_global = A_pos_global;
            if (L_straight > epsilon) {
                 node_pos_global += AF_vec_global * (s_arc_from_anchor / L_straight);
            }
            if (N_nodes_param[k_node_idx]) {
                N_nodes_param[k_node_idx]->XCurr.pos = node_pos_global;
                N_nodes_param[k_node_idx]->VCurr.v = Vec3(0., 0., 0.);
            }
        }
        return;
    }

    // ここからカテナリー理論による初期位置の決定
    double H_tension = 0.0;
    double p0_angle_out_val = 0.0;

    if (h_abs < epsilon) {
        if (L_APFP < epsilon) {
            silent_cerr("ModuleCatenaryLM::SetInitialValue: WARNING: Anchor and Fairlead coincide, L_orig > 0. Placing all nodes at anchor." << std::endl);
            for (unsigned int i = 0; i < Seg_param; ++i) {
                if (N_nodes_param[i]) {
                    N_nodes_param[i]->XCurr.pos = A_pos_global;
                    N_nodes_param[i]->VCurr.v = Vec3(0., 0., 0.);
                }
            }
            return;
        }

        double current_a = L_APFP;
        for (int iter = 0; iter < 50; ++iter) {
            if (current_a < epsilon) {
                current_a = epsilon;
            }
            double x_div_2a = L_APFP / (2.9 * current_a);
            double sinh_val = std::sinh(x_div_2a);
            double cosh_val = std::cosh(x_div_2a);
            double f_a = 2.0 * current_a * sinh_val - L_orig;
            double df_a = 2.0 * sinh_val - (L_APFP / current_a) * cosh_val;

            if (std::fabs(f_a) < epsilon*L_orig) break;
            if (std::fabs(df_a) < epsilon) break;

            current_a -= f_a / df_a;
            if (current_a <= 0) current_a = epsilon;
        }

        H_tension = current_a * w_orig;
        goto horizontal_catenary_solved;
    }

    {
        double d_geom = (L_APFP - (L_orig - h_abs)) / h_abs;
        double l_geom = L_orig / h_abs;
        double x_param_calc = 0.0;

        if (d_geom <= 0) {
            H_tension = 0.0;
            p0_angle_out_val = 0.0;
        } else {
            try {
                x_param_calc = ModuleCatenaryLM::rtsafe_catenary_local(epsilon, 1.e+8, xacc_orig, d_geom, l_geom, p0_angle_out_val);
                H_tension = x_param_calc * w_orig * h_abs;
            } catch (const ErrInterrupted& e) {
                silent_cerr("ModuleCatenaryLM::SetInitialValue: rtsafe_catenary_local failed: " << e.what() 
                            << " for d_geom=" << d_geom << ", l_geom=" << l_geom << std::endl);
                throw ErrGeneric("ModuleCatenaryLM::SetInitialValue: rtsafe solver failed.");
            }
        }
    }

    horizontal_catenary_solved:;
        std::vector<Vec3> local_corrds_nodes(Seg_param);
        if (std::fabs(H_tension) < epsilon && ! ((L_APFP - (L_orig - h_abs)) / h_abs <= 0 && h_abs > epsilon)) {
            silent_cerr("ModuleCatenaryLM::SetInitialValue: WARNING: H_tension is effectively zero. Assuming vertical drop / horizontal lay." << std::endl);
            H_tension = 0.0;
        }

        if (std::fabs(H_tension) < epsilon) { // H_tension = 0 (e.g. d_geom <= 0, or forced above)
            double L_on_seabed_or_slack_horizontal = 0.0;
            if (L_orig >= h_abs) {
                L_on_seabed_or_slack_horizontal = L_orig - h_abs;
            } else { // L_orig < h_abs, line can't even span vertically. This was checked by L_orig < L_straight.
                L_on_seabed_or_slack_horizontal = 0.0; 
            }
        
            for (unsigned int k_node_idx = 0; k_node_idx < Seg_param; ++k_node_idx) {
                double s_arc_from_anchor = L_orig * (static_cast<double>(Seg_param - k_node_idx) / static_cast<double>(Seg_param));
                double x_lp, z_lp;

                if (s_arc_from_anchor <= L_on_seabed_or_slack_horizontal + epsilon) { 
                    x_lp = s_arc_from_anchor;
                    z_lp = 0.0; 
                } else { // Node on vertical part
                    x_lp = L_on_seabed_or_slack_horizontal;
                    z_lp = s_arc_from_anchor - L_on_seabed_or_slack_horizontal;
                }
                local_coords_nodes[k_node_idx] = Vec3(x_lp, (h_FP_rel_AP >= 0 ? z_lp : -z_lp) , 0.0);
            }
        } else {
            double a = H_tension / w_orig;
            if (std::fabs(a) < epsilon) {
                silent_cerr("ModuleCatenaryLM::SetInitialValue: ERROR: Catenary parameter 'a' is near zero despite H_tension > 0. Check H_tension and w_orig." << std::endl);
                throw ErrGeneric("ModuleCatenaryLM::SetInitialValue: Catenary parameter 'a' is effectively zero.");
            }

            double x_F_lp = L_APFP;
            double z_F_lp = h_FP_rel_AP;

            double sigma_arg_den = 2.0 * a * std::sinh(x_F_lp / (2.0 * a));
            double cosh_sigma_arg_val;

            if (std::fabs(sigma_arg_den) < epsilon) {
                if (std::fabs(L_APFP) < epsilon) { // Indeed vertical line
                    if (std::fabs(L_orig - h_abs) > epsilon) {
                        silent_cerr("ModuleCatenaryLM::SetInitialValue: WARNING: Vertical line (L_APFP=0) with L_orig != h_abs. L_orig=" << L_orig << ", h_abs=" << h_abs <<". Results may be approximate for slack." << std::endl);
                    }
                    for (unsigned int k_node_idx = 0; k_node_idx < Seg_param; ++k_node_idx) {
                        double s_arc_from_anchor = L_orig * (static_cast<double>(Seg_param - k_node_idx) / static_cast<double>(Seg_param));
                        local_coords_nodes[k_node_idx] = Vec3(0.0, (h_FP_rel_AP >= 0 ? s_arc_from_anchor : -s_arc_from_anchor), 0.0);
                    }
                    goto transform_to_3d_nodes; 
                } else {
                    cosh_sigma_arg_val = 1.0; // Effectively straight line
                }
            } else {
                cosh_sigma_arg_val = L_orig / sigma_arg_den;
            }
        }

        double sigma = ModuleCatenaryLM::myacosh_local(local_robust_acosh_arg(cosh_sigma_arg_val));

        double avg_u;
        if (std::fabs(L_orig) < epsilon) avg_u = 0; // Avoid division by zero if L_orig is zero
        else avg_u = ModuleCatenaryLM::myatanh_local( std::max(-0.9999999, std::min(0.9999999, z_F_lp/L_orig)) );
        
        double u1 = avg_u - x_F_lp/(2.0*a);
        double x_v_lp = -a * u1; 
        double z_offset_lp = -a * ModuleCatenaryLM::mycosh_local(u1); 

        for (unsigned int k_node_idx = 0; k_node_idx < Seg_param; ++k_node_idx) {
            double s_arc_from_anchor = L_orig * (static_cast<double>(Seg_param - k_node_idx) / static_cast<double>(Seg_param));
            
            double sinh_U_val = s_arc_from_anchor / a + std::sinh(u1); 
            double U_target = ModuleCatenaryLM::myasinh_local(sinh_U_val);

            double x_lp = a * U_target + x_v_lp;
            double z_lp = a * ModuleCatenaryLM::mycosh_local(U_target) + z_offset_lp;
            local_coords_nodes[k_node_idx] = Vec3(x_lp, z_lp, 0.0);
        }
    
    transform_to_3d_nodes:;

        Vec3 u_horiz_global(0.0,0.0,0.0);
        if (L_APFP > epsilon) {
            u_horiz_global.dSet(1, AF_vec_global.dGet(1) / L_APFP);
            u_horiz_global.dSet(2, AF_vec_global.dGet(2) / L_APFP);
        } 
        
        for (unsigned int k_node_idx = 0; k_node_idx < Seg_param; ++k_node_idx) {
            if (N_nodes_param[k_node_idx]) {
                const Vec3& node_lp = local_coords_nodes[k_node_idx];
            
                Vec3 P_3D_node_global = A_pos_global;
                P_3D_node_global.dAdd(1, node_lp.dGet(1) * u_horiz_global.dGet(1)); 
                P_3D_node_global.dAdd(2, node_lp.dGet(1) * u_horiz_global.dGet(2)); 
                P_3D_node_global.dAdd(3, node_lp.dGet(2)); 

                N_nodes_param[k_node_idx]->XCurr.pos = P_3D_node_global;
                N_nodes_param[k_node_idx]->VCurr.v = Vec3(0., 0., 0.); // Zero initial velocity
            
            }
        }
        if (N_nodes_param[0]) {
            Vec3 calculated_F_pos = N_nodes_param[0]->GetXCurr().pos;
            Vec3 diff_F = calculated_F_pos - F_pos_global;
            if (diff_F.Norm() > epsilon * L_orig && L_orig > epsilon) {
                 silent_cerr("ModuleCatenaryLM::SetInitialValue: WARNING: Calculated fairlead position does not exactly match input. Diff norm: " << diff_F.Norm() << std::endl);
            }
        }
}


// Outputメソッド: シミュレーション結果の出力
void ModuleCatenaryLM::Output(OutputHandler& OH) const {
    // この要素が出力対象かどうかの判定
    if (bToBeOutput()) {
        // テキスト形式で出力するかどうかの判定
        if (OH.UseText(OutputHandler::LOADABLE)) {
            // フェアリーダーノードが存在し，有効かどうかの確認
            if (!N_nodes_param.empty() && N_nodes_param[0] != 0) {
                OH.Loadable() << GetLabel() // 要素のラベル
                              << " FairleadPos "
                              << N_nodes_param[0]->GetXCurr().dGet(1) << " " // フェアリーダーの X 座標
                              << N_nodes_param[0]->GetXCurr().dGet(2) << " " // フェアリーダーの Y 座標
                              << N_nodes_param[0]->GetXCurr().dGet(3) // // フェアリーダーの Z 座標
                              << std::endl; 
            } else {
                OH.Loadable() << GetLabel() << " Error: Fairlead node not available for output." << std::endl;
            }
            // ========== 各セグメントの張力や内部ノードの位置なども出力するならここに ============
        }
    }
}

// ================ ヤコビアン行列と残差ベクトルの次元を設定 ====================
// origin ではフェアリーダーポイントのみで 6 自由度で扱っていたが，ランプドマス法で内部ノードも全て管理する場合は大きく変わる
void ModuleCatenaryLM::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
    unsigned int num_nodes = Seg_param;
    unsigned int dof_per_node = 3;

    *piNumRows = num_nodes * dof_per_node;
    *piNumCols = num_nodes * dof_per_node;
}

void ModuleCatenaryLM::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
}

// ================ 残差 ==========================
SubVectorHandler& ModuleCatenaryLM::AssRes(
    SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr
) {

    integer iNumRows = 0;
    integer iNumCols = 0;
    WorkSpaceDim(&iNumRows, &iNumCols);

    WorkVec.ResizeReset(iNumRows);

    for (unsigned int i = 0; i < N_nodes_param.size(); ++i) {
        StructDispNode* node_i = N_nodes_param[i];
        if (node_i == 0) {
            continue;
        }

        const Vec3& pos_i = node_i->GetXCurr();
        const Vec3& vel_i = node_i->GetVCurr();

        unsigned int F_idx_start = i * 3;
        Vec3 total_force_on_node_i(0.0, 0.0, 0.0);

        // ====== 重力 ======
        doublereal node_mass_contribution = 0.0;
        if (i > 0) {
            if ((i-1) < P_param.size()) {
                node_mass_contribution += P_param[i-1].M_seg / 2.0;
            }
        }
        if (i < P_param.size()) {
            node_mass_contribution += P_param[i].M_seg / 2.0;
        }

        Vec3 gravity_force(0.0, 0.0, -node_mass_contribution * g_gravity_param);
        total_force_on_node_i += gravity_force;

        // ====== 弾性力 ======


        // ====== 減衰力 ======
        // 左側のセグメント (ノード i-1 と ノード i を接続) からの寄与
        if (i > 0) {
            StructDispNode* node_prev = N_nodes_param[i-1];
            if (node_prev && (i-1) < P_param.size() && P_param[i-1].CA_seg > 0.0) {
                const SegmentProperty& seg_prop_prev = P_param[i-1];
                const Vec3& pos_prev = node_prev->GetXCurr();
                const Vec3& vel_prev = node_prev->GetVCurr();

                Vec3 vec_prev_i = pos_i - pos_prev;
                doublereal length_prev_i = vec_prev_i.Norm();
                Vec3 unit_vec_prev_i(0.0, 0.0, 0.0);
                if (length_prev_i > 1.0e-12) {
                    unit_vec_prev_i = vec_prev_i / length_prev_i;
                }

                Vec3 rel_vel_prev_i = vel_i - vel_prev;
                doublereal axial_rel_vel = rel_vel_prev_i * unit_vec_prev_i;
                doublereal damping_force_magnitude = seg_prop_prev.CA_seg * axial_rel_vel;

                Vec3 damping_force_on_node_i_from_prev = unit_vec_prev_i * (-damping_force_magnitude);
                total_force_on_node_i += damping_force_on_node_i_from_prev;
            }
        }

        // 右側のセグメント (ノード i と ノード i+1 またはアンカー を接続) からの寄与
        if (i < P_param.size()) { // セグメントP_param[i]が存在する
            const SegmentProperty& seg_prop_curr = P_param[i];
            if (seg_prop_curr.CA_seg > 0.0) {
                Vec3 pos_partner(0.0,0.0,0.0);
                Vec3 vel_partner(0.0,0.0,0.0);
                bool partner_is_anchor = false;

                if (i < N_nodes_param.size() - 1) {
                    StructDispNode* node_next = N_nodes_param[i+1];
                    if (node_next) {
                        pos_partner = node_next->GetXCurr();
                        vel_partner = node_next->GetVCurr();
                    } else {
                        continue;
                    }
                } else {
                    pos_partner = Vec3(APx_orig, APy_orig, APz_orig);
                    vel_partner = Vec3(0.0, 0.0, 0.0);
                    partner_is_anchor = true;
                }

                Vec3 vec_i_partner = pos_partner - pos_i;
                doublereal length_i_partner = vec_i_partner.Norm();
                Vec3 unit_vec_i_partner(0.0, 0.0, 0.0);
                if (length_i_partner > 1.0e-12) {
                    unit_vec_i_partner = vec_i_partner / length_i_partner;
                }

                Vec3 rel_vel_i_partner = vel_partner - vel_i;
                doublereal axial_rel_vel = rel_vel_i_partner * unit_vec_i_partner;
                doublereal damping_force_magnitude = seg_prop_curr.CA_seg * axial_rel_vel;

                Vec3 damping_force_on_node_i_from_curr_seg =  unit_vec_i_partner * (-damping_force_magnitude);
                total_force_on_node_i += damping_force_on_node_i_from_curr_seg;
            }
        }

        // ====== WorkVec への格納 ======
        WorkVec.PutCoef(F_idx_start + 1, total_force_on_node_i.dGet(1));
        WorkVec.PutCoef(F_idx_start + 2, total_force_on_node_i.dGet(2));
        WorkVec.PutCoef(F_idx_start + 3, total_force_on_node_i.dGet(3));
    } 

    return WorkVec;
}

// 初期形状をカテナリー理論で決める
SubVectorHandler& ModuleCatenaryLM::InitialAssRes(
    SubVectorHandler& WorkVec,
    const VectorHandler& XCurr
) {
    return WorkVec;
}

// ============== 全体ヤコビアン ====================
VariableSubMatrixHandler& ModuleCatenaryLM::AssJac(
    VariableSubMatrixHandler& WorkMat,
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr
) {
    return WorkMat;
}

VariableSubMatrixHandler& ModuleCatenaryLM::InitialAssJac(
    VariableSubMatrixHandler& WorkMat,
    const VectorHandler& /*X*/
) {
    WorkMat.SetNullMatrix();
    return WorkMat;
}

// この要素がMBDynのモデル内で接続している，あるいは管理しているノードの数を返す
// MBDyn は要素とノード間の接続情報を構築・管理するためにこの情報を利用する
unsigned int ModuleCatenaryLM::iGetNumConnectedNodes(void) const {
    return Seg_param; // これで完成のはず
}

// 後回し
std::ostream& ModuleCatenaryLM::Restart(std::ostream& out) const {
    out << "# ModuleCatenaryLM (Label: " << GetLabel() << ") Restart: Not implemented yet." << std::endl;
    return out;
}

// 初期化フェーズで特別な解析をするなら，初期のみ自由度を付与する
unsigned int ModuleCatenaryLM::iGetInitialNumDof(void) const {
    return 0;
}


void ModuleCatenaryLM::SetValue(
    DataManager* /*pDM*/,
    VectorHandler& /*X*/,
    VectorHandler& /*XP*/,
    SimulationEntity::Hints* /*pHints*/ // 線形化が必要かどうか
)
{
    // 重力荷重は AssRes 内で扱うか，要素外で扱うのでここは空実装
}

// ==== ModuleCatenaryLM 要素が管理している i 番目のノードへのポインタを返すためのインターフェース ====
// MBDyn のソルバーや他のモジュールが，この要素に関連付けられたノードの情報（位置・速度・力など）にアクセスしたり，ノードの状態を変更したりするために使用する
// const 版：この関数自身が const であり，この関数を呼び出しても ModuleCatenaryLM オブジェクトの状態は変更されない
const Node* ModuleCatenaryLM::pGetNode(unsigned int i) const {
    if (i < N_nodes_param.size() && N_nodes_param[i] != 0) {
        return N_nodes_param[i];
    }
    return 0;
}

// 非 const 版：ノードの状態を変更することが可能
Node* ModuleCatenaryLM::pGetNode(unsigned int i) {
    if (i < N_nodes_param.size() && N_nodes_param[i] != 0) {
        return N_nodes_param[i];
    }
    return 0;
}

bool catenary_lm_set(void) {
#ifdef DEBUG
    std::cerr << __FILE__ << ":" << __LINE__ << ": bool catenary_lm_set(void)" << std::endl;
#endif
    UserDefinedElemRead *rf = new UDERead<ModuleCatenaryLM>; // 新しいクラス名でテンプレート特殊化
    if (!SetUDE("catenary_lm", rf)) {
        delete rf;
        return false;
    }
    return true;
}

#ifndef STATIC_MODULES // MBDyn の標準的な動的モジュールロードの仕組み
extern "C" {
    int module_init(const char *module_name, void *pdm /*DataManager* */, void *php /* MBDynParser* */) {
        if (!catenary_lm_set()) { // 新しいセット関数を呼ぶ
            return -1; // 失敗
        }
        return 0; // 成功
    }
} // extern "C"
#endif // ! STATIC_MODULES
