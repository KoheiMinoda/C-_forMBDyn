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
// #include "vec3.h"          // Vec3 (座標や力のベクトル表現に必須)
// #include "mbdynparser.h"   // MBDynParser (パラメータ読み込みに必須)
// #include "output_handler.h"// OutputHandler (結果出力に必須)
#include "node.h"          // Node (pGetNode の戻り値の基底クラス)

// 追加の標準ライブラリ
#include <vector>       // std::vector (ノードやセグメントのリスト管理用)

class ModuleCatenaryLM : virtual public Elem, public UserDefinedElem
{
public:
    ModuleCatenaryLM(unsigned uLabel, const DofOwner *pD0, DataManager* pDM, MBDynParser& HP);
    virtual ~ModuleCatenaryLM(void);

    virtual void Output(OutputHandler& OH) const;
    virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

    VariableSubMatrixHandler&
    AssJac(VariableSubMatrixHandler& WorkMat, 
        doublereal dCoef, 
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr
    );

    SubVectorHandler&
    AssRes(SubVectorHandler& WorkVec, 
        doublereal dCoef, 
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr
    );

    virtual unsigned int iGetNumConnectedNodes(void) const;
    virtual void SetInitialValue(VectorHandler& X, VectorHandler& XP);
    virtual std::ostream& Restart(std::ostream& out) const;

    virtual unsigned int iGetInitialNumDof(void) const;
    virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    virtual VariableSubMatrixHandler&
    InitialAssJac(VariableSubMatrixHandler& WorkMat, const VectorHandler& XCurr);
    virtual SubVectorHandler&
    InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

    virtual const Node* pGetNode(unsigned int i) const;
    virtual Node* pGetNode(unsigned int i);

private:
    double APx_orig; // アンカーの x 座標
    double APy_orig; // アンカーの y 座標
    double APz_orig; // アンカーの z 座標
    double L_orig; // ライン全長
    double w_orig; // チェーンの単重
    double xacc_orig; // rtsafe 関数で使用するパラメータ
    DriveOwner FSF_orig; // 元の Force Scale Factor

    static double myasinh_local(double val);
    static double myacosh_local(double val);
    static double myatanh_local(double val);
    static void funcd_catenary_local(double x_param, double xacc, double &f_val, double &df_val, double d_geom, double l_geom, double &p0_angle);
    static double rtsafe_catenary_local(double x1_bounds, double x2_bounds, double xacc_tol, double d_geom, double l_geom, double &p0_angle_out);

    unsigned int Seg_param;

    // ノード配列 (StructNode* の動的配列)
    // N_nodes_param[0] はフェアリーダーノード
    // N_nodes_param[1] ... N_nodes_param[Seg_param-1] は内部質量点ノード
    // アンカーは固定座標(APx_orig, APy_orig, APz_orig)として別途扱う
    // このベクトルサイズは Seg_param になる
    std::vector<StructNode*> N_nodes_param;

    // 各セグメントの特性
    struct SegmentProperty {
        doublereal L0_seg;  // セグメントの初期自然長
        doublereal M_seg;   // この「セグメント」の質量
        doublereal EA_seg;  // セグメントの軸剛性 (EA)
        doublereal CA_seg;  // セグメントの軸方向減衰係数 (オプション)

        SegmentProperty() : L0_seg(0.0), M_seg(0.0), EA_seg(0.0), CA_seg(0.0) {}
    };
    std::vector<SegmentProperty> P_param; // サイズは Seg_param (セグメント数)

    doublereal g_gravity_param;
};

ModuleCatenaryLM::ModuleCatenaryLM(
    unsigned uLabel,
    const DofOwner *pD0,
    DataManager* pDM,
    MBDynParser& HP
)
    : Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pD0),
        APx_orig(0.0), APy_orig(0.0), APz_orig(0.0),
        L_orig(0.0), w_orig(0.0), xacc_orig(1e-6),
        Seg_param(20),
        g_gravity_param(9.80665)
    {
        if (HP.IsKeyWord("help")) {
            silent_cout(
                "\n"
                "Module:     ModuleCatenaryLM (Lumped Mass Catenary Mooring Line - Segments Fixed to 20)\n"
                "Usage:      catenary_lm, label,\n"
                "                fairlead_node_label,\n" // アンカーは固定座標、セグメント数は固定
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

        // --- パラメータ読み込み開始 ---
        // 読み込み順は元の module-catenary_origin.cc とヘルプメッセージに合わせる

        // 1. フェアリーダーノードのラベル (StructNode*)
        unsigned int fairlead_node_label = HP.GetInt();

        // 2. 全長 L_orig
        L_orig = HP.GetReal();
        if (L_orig <= 0.0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        // 3. 単位重量 w_orig
        w_orig = HP.GetReal();
        if (w_orig <= 0.0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        // 4. rtsafeの計算精度 xacc_orig
        xacc_orig = HP.GetReal();
        if (xacc_orig <= 0.0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        // 5. アンカー固定座標 APx_orig, APy_orig, APz_orig
        APx_orig = HP.GetReal();
        APy_orig = HP.GetReal();
        APz_orig = HP.GetReal();

        // 6. 軸剛性 EA (ランプドマス法用)
        doublereal EA_val;
        if (HP.IsKeyWord("EA")) { // キーワードがある方がより頑健
            EA_val = HP.GetReal();
        } else {
            EA_val = HP.GetReal(); // キーワードなしで順番に読む場合
        }
        if (EA_val <= 0.0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        // 7. 軸方向減衰係数 CA (オプション)
        doublereal CA_val = 0.0; // デフォルトは減衰なし
        if (HP.IsKeyWord("CA")) {
            CA_val = HP.GetReal();
            if (CA_val < 0.0) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        // 8. 重力加速度 g (オプション)
        if (HP.IsKeyWord("gravity")) {
            g_gravity_param = HP.GetReal();
            if (g_gravity_param < 0.0) { // 通常は正の値
                // throw ErrGeneric(MBDYN_EXCEPT_ARGS); // エラーにするか、警告にとどめるか
            }
        }

        // 9. FSFの読み込み (元のコードから流用)
        // 元のコードでは "Force" "scale" "factor" とスペース区切りでキーワードを認識しようとしている。
        if (HP.IsKeyWord("force_scale_factor")) { // より一般的なキーワード形式
            FSF_orig.Set(HP.GetDriveCaller());
        } else if (HP.IsKeyWord("Force") && HP.IsKeyWord("scale") && HP.IsKeyWord("factor")) { // 元の形式を試みる
            FSF_orig.Set(HP.GetDriveCaller());
        } else {
            FSF_orig.Set(new OneDriveCaller); // デフォルトはスケール1.0
        }

        // 10. 出力フラグの設定 (元のコードから流用)
        SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

        // --- ノードの確保と初期化 ---
        // N_nodes_param にはフェアリーダーと Seg_param-1 個の内部ノードを格納
        // 合計 Seg_param 個の StructNode* を持つ
        N_nodes_param.resize(Seg_param);

        // フェアリーダーノードの取得と設定
        N_nodes_param[0] = dynamic_cast<StructNode*>(pDM->pFindNode(fairlead_node_label));
        if (N_nodes_param[0] == 0) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
        AddNode(N_nodes_param[0]);

        Vec3 temp_initial_pos(0.0, 0.0, 0.0);

        for (unsigned int i=1; i<Seg_param; ++i) {
            N_nodes_param[i] = pDM->CreateStructNode(pD0, temp_initial_pos, Vec3(0.0, 0.0, 0.0)); // 初期速度 0
            if (N_nodes_param[i] == 0) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
            AddNode(N_nodes_param[i]);
        }

        P_param.resize(Seg_param); // Seg_param 個のセグメント
        doublereal L0_s = L_orig / static_cast<doublereal>(Seg_param); // 各セグメントの自然長
        doublereal mass_per_segment = (w_orig / g_gravity_param) * L0_s; // 各セグメントの質量

        for (unsigned int i = 0; i < Seg_param; ++i) {
            P_param[i].L0_seg = L0_s;
            P_param[i].M_seg = mass_per_segment;
            P_param[i].EA_seg = EA_val;
            P_param[i].CA_seg = CA_val;
        }

        // ログ出力
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

void ModuleCatenaryLM::funcd_catenary_local(double x_param, double xacc, double &f_val, double &df_val, double d_geom, double l_geom, double &p0_angle) {
    int i_internal_iter,max_internal_iter;
    double f1_internal, df1_internal;
    max_internal_iter = 1000;

    if(x_param==0.0) {
        f_val=-d_geom;
        df_val=0.0;
        p0_angle=0.0;
    } else if(x_param>0.0) {
        if(l_geom<=0.0){
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
        } else {
            if(x_param>(l_geom*l_geom-1.0)/2.0) {
                p0_angle=0.0;
                for(i_internal_iter=0; i_internal_iter<max_internal_iter; i_internal_iter++) {
                    double cos_p0 = std::cos(p0_angle);
                    if (std::fabs(cos_p0) < 1e-9) { df1_internal = 1.0; break; } // 保護
                    double func1_internal = 1.0/x_param+1.0/cos_p0;
                    double term_in_sqrt_f1 = func1_internal*func1_internal-1.0;
                    if (term_in_sqrt_f1 < 0.0) term_in_sqrt_f1 = 0.0;

                    f1_internal=x_param*(std::sqrt(term_in_sqrt_f1)-std::tan(p0_angle))-l_geom;

                    if (std::fabs(cos_p0) < 1e-9 || term_in_sqrt_f1 < 1e-12 ) { df1_internal = 1.0; break; } // 保護
                    df1_internal = x_param * (func1_internal * std::tan(p0_angle) / (cos_p0 * std::sqrt(term_in_sqrt_f1)) - (std::tan(p0_angle)*std::tan(p0_angle)) - 1.0);

                    if (std::fabs(df1_internal) < 1e-9) { break; }
                    p0_angle=p0_angle-f1_internal/df1_internal;

                    cos_p0 = std::cos(p0_angle);
                    if (std::fabs(cos_p0) < 1e-9) { break; }
                    func1_internal = 1.0/x_param+1.0/cos_p0;
                    term_in_sqrt_f1 = func1_internal*func1_internal-1.0;
                    if (term_in_sqrt_f1 < 0.0) term_in_sqrt_f1 = 0.0;
                    f1_internal=x_param*(std::sqrt(term_in_sqrt_f1)-std::tan(p0_angle))-l_geom;

                    if(std::fabs(f1_internal)<xacc) { break; }
                }
                if(i_internal_iter==max_internal_iter && std::fabs(f1_internal)>xacc) {
                }

                double X_2_internal = l_geom/x_param+std::tan(p0_angle);
                double X_3_internal = std::tan(p0_angle);
                f_val=x_param*(myasinh_local(X_2_internal)-myasinh_local(X_3_internal))-l_geom+1.0-d_geom;
                
                double term_in_sqrt_df = X_2_internal*X_2_internal+1.0;
                // if (term_in_sqrt_df < 0.0) term_in_sqrt_df = 0.0; // asinhの引数は実数なので常に正
                if (std::fabs(x_param * std::sqrt(term_in_sqrt_df)) < 1e-12 ) { df_val = 1.0; } // 保護
                else {
                    df_val=myasinh_local(X_2_internal)-myasinh_local(X_3_internal)-l_geom/(x_param*std::sqrt(term_in_sqrt_df));
                }
            } else {
                double X_5_internal = 1.0/x_param+1.0;
                if (X_5_internal < 1.0) X_5_internal = 1.0; // acoshの引数保護
                double term_sqrt1 = 1.0+2.0*x_param;
                if (term_sqrt1 < 0.0) term_sqrt1 = 0.0;
                double term_sqrt2 = X_5_internal*X_5_internal-1.0;
                if (term_sqrt2 < 0.0) term_sqrt2 = 0.0;

                f_val=x_param*myacosh_local(X_5_internal)-std::sqrt(term_sqrt1)+1.0-d_geom;
                if (std::fabs(x_param * std::sqrt(term_sqrt2)) < 1e-12 || std::fabs(std::sqrt(term_sqrt1)) < 1e-12) { // 保護
                     df_val = myacosh_local(X_5_internal); // 近似
                } else {
                     df_val=myacosh_local(X_5_internal)-1.0/std::sqrt(term_sqrt1)-1.0/(x_param*std::sqrt(term_sqrt2));
                }
                p0_angle=0.0;
            }
        }
    } else {
        throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
    }
}

double ModuleCatenaryLM::rtsafe_catenary_local(double x1_bounds, double x2_bounds, double xacc_tol, double d_geom, double l_geom, double &p0_angle_out) {
    const int MAXIT_internal=1000;
    int j_internal;
    double fh_internal,fl_internal,xh_internal,xl_internal;
    double dx_internal,dxold_internal,f_internal,temp_internal,rts_internal;
    double p1_internal, p2_internal;
    double df_not_used, df_internal; // funcd_catenary_local がdfを要求するため

    funcd_catenary_local(x1_bounds, xacc_tol, fl_internal, df_not_used, d_geom, l_geom, p1_internal);
    funcd_catenary_local(x2_bounds, xacc_tol, fh_internal, df_not_used, d_geom, l_geom, p2_internal);

    if((fl_internal>0.0 && fh_internal>0.0) || (fl_internal<0.0 && fh_internal<0.0)) {
        throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
    }
    if(fl_internal==0.0) {
        p0_angle_out  = p1_internal;
        return x1_bounds;
    }
    if(fh_internal==0.0) {
        p0_angle_out  = p2_internal;
        return x2_bounds;
    }

    if(fl_internal<0.0) {
        xl_internal=x1_bounds;
        xh_internal=x2_bounds;
    } else {
        xh_internal=x1_bounds;
        xl_internal=x2_bounds;
    }

    rts_internal=0.5*(x1_bounds+x2_bounds);
    dxold_internal=std::fabs(x2_bounds-x1_bounds);
    dx_internal=dxold_internal;
    funcd_catenary_local(rts_internal, xacc_tol, f_internal, df_internal, d_geom, l_geom, p0_angle_out);

    for(j_internal=0; j_internal<MAXIT_internal; j_internal++) {
        if((((rts_internal-xh_internal)*df_internal-f_internal)*((rts_internal-xl_internal)*df_internal-f_internal)>0.0)
           || (std::fabs(2.0*f_internal) > std::fabs(dxold_internal*df_internal))) {
            dxold_internal  = dx_internal;
            dx_internal     =0.5*(xh_internal-xl_internal);
            rts_internal    =xl_internal+dx_internal;
            if(xl_internal==rts_internal){ return rts_internal; }
        } else {
            dxold_internal=dx_internal;
            if (std::fabs(df_internal) < 1e-12) { // ゼロ除算を避ける
                dx_internal     =0.5*(xh_internal-xl_internal); // 二分法にフォールバック
                rts_internal    =xl_internal+dx_internal;
                if(xl_internal==rts_internal){ return rts_internal; }
            } else {
                dx_internal=f_internal/df_internal;
            }
            temp_internal=rts_internal;
            rts_internal-=dx_internal;
            if(temp_internal==rts_internal) {return rts_internal;}
        }

        if(std::fabs(dx_internal)<xacc_tol) { return rts_internal; }

        funcd_catenary_local(rts_internal, xacc_tol, f_internal, df_internal, d_geom, l_geom, p0_angle_out);

        if(f_internal<0.0){
            xl_internal=rts_internal;
        } else {
            xh_internal=rts_internal;
        }
    }

    throw ErrInterrupted(MBDYN_EXCEPT_ARGS);
}


void ModuleCatenaryLM::SetInitialValue(VectorHandler& X, VectorHandler& XP) {

    if (GetLabel() != 0) { // GetLabel() が有効な値を持つことを期待 (通常はMBDynが設定)
        std::cout << "ModuleCatenaryLM (Label: " << GetLabel()
                  << "): SetInitialValue called. Initializing node positions with linear interpolation."
                  << std::endl;
    } else {
        std::cout << "ModuleCatenaryLM (Label: N/A): SetInitialValue called. Initializing node positions with linear interpolation."
                  << std::endl;
    }


    // 1. フェアリーダーノードの初期位置を取得
    if (N_nodes_param.empty() || N_nodes_param[0] == 0) {
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    // SetInitialValue が呼ばれる時点では、接続ノードの初期位置は MBDyn によって
    // グローバル状態ベクトル X に設定されていると期待されます。
    Vec3 fairlead_pos = N_nodes_param[0]->GetXCurr();

    // アンカーの固定座標
    Vec3 anchor_pos(APx_orig, APy_orig, APz_orig);

    // --- 各ランプドマス点の初期位置を直線補間で設定 ---

    // N_nodes_param[0] (フェアリーダー) の初期位置と速度を設定
    // (GetXCurr(X) で取得した値で再度設定する形になるが、速度を確実に0にする意味もある)
    N_nodes_param[0]->PutXCurr(fairlead_pos, X);
    N_nodes_param[0]->PutVCurr(Vec3(0.0, 0.0, 0.0), XP); // 初期速度ゼロ

    if (Seg_param < 1) { // このチェックはコンストラクタでも行っているが念のため
        return; // 何もせず終了
    }

    for (unsigned int i = 1; i < Seg_param; ++i) { // N_nodes_param のインデックスに対応
                                                  // i は内部ノードのインデックス (1 から Seg_param-1)
        // ノード i は、フェアリーダーから見て (i / Seg_param) の割合の位置にある
        doublereal ratio = static_cast<doublereal>(i) / static_cast<doublereal>(Seg_param);
        Vec3 node_pos = fairlead_pos + (anchor_pos - fairlead_pos) * ratio;

        if (N_nodes_param[i] == 0) { // 安全チェック
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
        N_nodes_param[i]->PutXCurr(node_pos, X);
        N_nodes_param[i]->PutVCurr(Vec3(0.0, 0.0, 0.0), XP); // 初期速度ゼロ

        // デバッグ用出力 (オプション)
        // std::cout << "  Node " << i << " (Internal) initial position: "
        //           << node_pos.dGet(1) << ", " << node_pos.dGet(2) << ", " << node_pos.dGet(3) << std::endl;
    }

    // アンカーポイントの処理：
    // アンカーは固定座標 (APx_orig, APy_orig, APz_orig) であり、MBDynの動的ノードではないため、
    // X や XP にアンカーの情報を設定する必要はありません。
    // N_nodes_param ベクターには、アンカーに対応する StructNode* は格納していません。

    // ログ出力 (オプション)
    std::cout << "ModuleCatenaryLM (Label: " << GetLabel()
              << "): Initial positions set for " << Seg_param << " nodes (including fairlead, excluding fixed anchor)."
              << std::endl;
    std::cout << "  Fairlead at: (" << fairlead_pos.dGet(1) << ", " << fairlead_pos.dGet(2) << ", " << fairlead_pos.dGet(3) << ")" << std::endl;
    if (Seg_param > 1 && N_nodes_param[Seg_param-1] != 0) { // 最後の内部ノード (アンカーの直前)
        Vec3 last_internal_node_pos = N_nodes_param[Seg_param-1]->GetXCurr(X);
         std::cout << "  Last internal node (" << Seg_param-1 << ") at: ("
                   << last_internal_node_pos.dGet(1) << ", " << last_internal_node_pos.dGet(2) << ", " << last_internal_node_pos.dGet(3) << ")" << std::endl;
    }
    std::cout << "  Anchor (fixed) at: (" << anchor_pos.dGet(1) << ", " << anchor_pos.dGet(2) << ", " << anchor_pos.dGet(3) << ")" << std::endl;

}

// (前回までのインクルード、クラス定義、コンストラクタ、静的ヘルパー関数、SetInitialValueの実装は変更なし)
// ...

// ---- ModuleCatenaryLM クラスの他のメソッドのスタブ実装 ----

// Outputメソッド: シミュレーション結果の出力
void ModuleCatenaryLM::Output(OutputHandler& OH) const {
    if (bToBeOutput()) { // MBDynの標準的な出力判定
        if (OH.UseText(OutputHandler::LOADABLE)) {
            // 簡単な例として、要素ラベルとフェアリーダーノードの現在の座標を出力
            if (!N_nodes_param.empty() && N_nodes_param[0] != 0) {
                OH.Loadable() << GetLabel() // 要素のラベル
                              << " FairleadPos "
                              << N_nodes_param[0]->GetXCurr().dGet(1) << " "
                              << N_nodes_param[0]->GetXCurr().dGet(2) << " "
                              << N_nodes_param[0]->GetXCurr().dGet(3)
                              << std::endl;
            } else {
                OH.Loadable() << GetLabel() << " Error: Fairlead node not available for output." << std::endl;
            }
            // 将来的には、各セグメントの張力や内部ノードの位置なども出力対象になります。
        }
        // 他の出力形式 (バイナリなど) のサポートもここに追加できます。
    }
}

// WorkSpaceDimメソッド: ヤコビアン行列と残差ベクトルの次元を設定
void ModuleCatenaryLM::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
    *piNumRows = 0;
    *piNumCols = 0;
}

// AssJacメソッド: ヤコビアン行列の計算 (スタブ)
VariableSubMatrixHandler& ModuleCatenaryLM::AssJac(
    VariableSubMatrixHandler& WorkMat,
    doublereal dCoef,
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr
) {
    if (WorkMat.pGetSubMatrix() != 0) { // 有効な行列ハンドラかチェック
        WorkMat.SetNullMatrix(); // または WorkMat.Zero(); (サイズが設定されていれば)
    }
    // silent_cout << "ModuleCatenaryLM (Label: " << GetLabel() << "): AssJac called (stub)." << std::endl;
    return WorkMat;
}

// AssResメソッド: 残差ベクトルの計算 (スタブ)
SubVectorHandler& ModuleCatenaryLM::AssRes(
    SubVectorHandler& WorkVec,
    doublereal dCoef,
    const VectorHandler& XCurr,
    const VectorHandler& XPrimeCurr
) {
    return WorkVec;
}

// iGetNumConnectedNodesメソッド: 接続されているノード数を返す
unsigned int ModuleCatenaryLM::iGetNumConnectedNodes(void) const {
    // N_nodes_param にはフェアリーダーと Seg_param-1 個の内部ノードが格納されている。
    // 合計 Seg_param 個のノードをこの要素が「接続」しているとMBDynに伝える。
    // これらは AddNode() で登録されたノード。
    return Seg_param;
}

// Restartメソッド: リスタート情報の出力 (スタブ)
std::ostream& ModuleCatenaryLM::Restart(std::ostream& out) const {
    out << "# ModuleCatenaryLM (Label: " << GetLabel() << ") Restart: Not implemented yet." << std::endl;
    // 将来的には、必要な内部状態（もしあれば）をここに出力します。
    // 通常、ノードの状態はMBDynが保存し、要素のパラメータは入力ファイルから再読み込みされます。
    return out;
}

// iGetInitialNumDofメソッド: 初期化時に使われる自由度数 (スタブ)
unsigned int ModuleCatenaryLM::iGetInitialNumDof(void) const {
    // SetInitialValue を使う場合、これが0でもMBDynが呼び出すことがあります。
    // 要素が独自の初期化自由度を持たない場合は0。
    return 0;
}

// InitialWorkSpaceDimメソッド: 初期化時のワークスペース次元 (スタブ)
void ModuleCatenaryLM::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
    // iGetInitialNumDof() が0なら、通常これも0。
    *piNumRows = 0;
    *piNumCols = 0;
}

// InitialAssJacメソッド: 初期化時のヤコビアン (スタブ)
VariableSubMatrixHandler& ModuleCatenaryLM::InitialAssJac(
    VariableSubMatrixHandler& WorkMat,
    const VectorHandler& XCurr
) {
    // InitialWorkSpaceDim が 0,0 なら、呼ばれないか、何もしない。
    if (WorkMat.pGetSubMatrix() != 0) {
        WorkMat.SetNullMatrix();
    }
    // ASSERT(0); // 元のコードにあったように、呼ばれるべきではない場合
    return WorkMat;
}

// InitialAssResメソッド: 初期化時の残差 (スタブ)
SubVectorHandler& ModuleCatenaryLM::InitialAssRes(
    SubVectorHandler& WorkVec,
    const VectorHandler& XCurr
) {
    // InitialWorkSpaceDim が 0,0 なら、呼ばれないか、何もしない。
    // if (WorkVec.iGetSize() > 0) { // サイズが設定されていればリセット
    //     WorkVec.ResizeReset(0);
    // }
    // ASSERT(0); // 元のコードにあったように、呼ばれるべきではない場合
    return WorkVec;
}

// pGetNodeメソッド (const版): 指定されたインデックスのノードへのポインタを返す
const Node* ModuleCatenaryLM::pGetNode(unsigned int i) const {
    if (i < N_nodes_param.size() && N_nodes_param[i] != 0) {
        return N_nodes_param[i];
    }
    return 0; // または適切なエラー処理
}

// pGetNodeメソッド (非const版): 指定されたインデックスのノードへのポインタを返す
Node* ModuleCatenaryLM::pGetNode(unsigned int i) {
    if (i < N_nodes_param.size() && N_nodes_param[i] != 0) {
        return N_nodes_param[i];
    }
    return 0; // または適切なエラー処理
}


// ---- モジュール登録関数の実装 (ファイル末尾に配置) ----
// (bool catenary_lm_set(void) と extern "C" int module_init(...) は前回提示した通り)
bool catenary_lm_set(void) {
#ifdef DEBUG
    // MBDynのデバッグ出力マクロがあればそれを使う方が良い
    std::cerr << __FILE__ << ":" << __LINE__ << ": bool catenary_lm_set(void)" << std::endl;
#endif
    UserDefinedElemRead *rf = new UDERead<ModuleCatenaryLM>; // 新しいクラス名でテンプレート特殊化

    // MBDynに登録する要素名を新しい名前に変更 (例: "catenary_lm")
    // 元のSetUDEはdataman.hで定義されているか、MBDynのグローバル関数のはず
    if (!SetUDE("catenary_lm", rf)) {
        delete rf;
        return false;
    }
    return true;
}

#ifndef STATIC_MODULES // MBDynの標準的な動的モジュールロードの仕組み
extern "C" {
    int module_init(const char *module_name, void *pdm /*DataManager* */, void *php /* MBDynParser* */) {
        if (!catenary_lm_set()) { // 新しいセット関数を呼ぶ
            return -1; // 失敗
        }
        return 0; // 成功
    }
} // extern "C"
#endif // ! STATIC_MODULES
